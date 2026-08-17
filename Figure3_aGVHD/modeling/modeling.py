import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, roc_curve

# Load Data
X = pd.read_csv("./modeling_counts_cpm_tmm_log.csv", index_col=0).T
metadata = pd.read_csv("./modeling_metadata.csv", index_col="sample_id")

# Extract y (expGroup) and ORIGIN
y = metadata["expGroup"]
origin = metadata["ORIGIN"]

# Create a stratification column combining expGroup & ORIGIN
metadata["stratify_group"] = metadata["expGroup"].astype(str) + "_" + metadata["ORIGIN"].astype(str)


# Store results
auc_results = []

# Run 100 iterations with different train-test splits
for seed in range(101):
# for seed in range(10):


    print(f"Iteration {seed+1}/100...")

    # Step 1: Stratified train-test split based on expGroup & ORIGIN
    cpm_train, cpm_test, y_train, y_test, strat_train, strat_test = train_test_split(
        X, y, metadata["stratify_group"], test_size=0.2, stratify=metadata["stratify_group"], random_state=seed
    )

    # Step 2: Remove features where at least 75% of samples in train set have values ≤ 0.5
    mask = (cpm_train > 0.5).sum(axis=0) >= (0.75 * len(cpm_train))  # Keep features where ≥ 25% of values are > 0.5
    cpm_train_filtered = cpm_train.loc[:, mask]

    # Step 3: Subset test set to selected features
    cpm_test_filtered = cpm_test.loc[:, cpm_train_filtered.columns]

    # Step 4: Define LASSO model and cross-validation
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=seed)
    lasso = LogisticRegression(penalty="l1", class_weight="balanced", max_iter=int(1e6), solver="liblinear", random_state=seed)
    lasso_cv = GridSearchCV(lasso, param_grid={"C": np.logspace(-2, 2, 30)}, scoring="roc_auc", cv=cv, n_jobs=-1)

    en = LogisticRegression(
        penalty='elasticnet',
        solver='saga',
        class_weight='balanced',
        max_iter=int(1e3),
        random_state=42
    )
    en_cv = GridSearchCV(en, param_grid={"C": np.logspace(-2, 1, 5), "l1_ratio": [.5, .7, .9]}, scoring="roc_auc", cv=cv, n_jobs=-1)

    # Step 5: Preprocessing + Model Pipeline
    preprocessing = Pipeline([
        ("variance", VarianceThreshold(0.5)),
        ("std", StandardScaler()),
        ("impute", SimpleImputer(strategy="median"))
    ])
    
    pipeline = Pipeline([
        # ("preprocessing", preprocessing),
        ("lasso_cv", lasso_cv)
    ])

    # Train the model
    pipeline.fit(cpm_train_filtered, y_train)

    # Step 6: Get best model and make predictions
    best_lasso = pipeline.named_steps["lasso_cv"].best_estimator_
    test_probs = best_lasso.predict_proba(cpm_test_filtered)[:, 1]

    # Calculate AUC and store results
    test_auc = roc_auc_score(y_test, test_probs)
    auc_results.append(test_auc)

    # train_probs = best_lasso.predict_proba(cpm_train_filtered)[:, 1]
    # train_auc = roc_auc_score(y_train, train_probs)
    # auc_results.append(train_auc)


# Convert results to DataFrame
auc_df = pd.DataFrame({"Iteration": np.arange(1, 102), "Test_AUC": auc_results})
# auc_df = pd.DataFrame({"Iteration": np.arange(1, 11), "Test_AUC": auc_results})


# Save results to CSV
auc_df.to_csv("lasso_test_auc_results.csv", index=False)

# Print summary statistics
print("\nTest AUC Summary:")
print(auc_df["Test_AUC"].describe())

# Plot AUC Distribution
plt.figure(figsize=(8, 6))
sns.histplot(auc_df["Test_AUC"], bins=20, kde=True)
plt.xlabel("Test AUC")
plt.ylabel("Frequency")
plt.title("Distribution of Test AUCs across 100 Iterations")
plt.show()