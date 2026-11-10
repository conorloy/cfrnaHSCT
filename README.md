# cfrnaHSCT

Analysis code for our study of cell-free RNA (cfRNA) in patients undergoing hematopoietic stem cell transplantation (HSCT), including trajectories after transplant and graft-versus-host disease (GVHD).

## Repository structure

| Folder | Contents |
| --- | --- |
| `1_DATA/` | Data preparation notebooks (data files not included) |
| `Figure1_projectOverview/` | Cohort and sample overview |
| `Figure2_trajectories/` | cfRNA trajectories after transplant |
| `Figure3_aGVHD/` | Acute GVHD analysis and modeling |
| `Figure4_LGVHD/` | Late GVHD analysis and modeling |
| `Figure5_MCCTaper/` | MCC taper cohort analysis and modeling |
| `Figure6_overlap/` | Overlap analysis across cohorts and tissue mapping |
| `Figure_Supplement/` | Extended Data figures |

Shared R helpers are in `functions.R` and `theme_ggplot_cfrna.R`. The full file tree is in `structure.txt`.
