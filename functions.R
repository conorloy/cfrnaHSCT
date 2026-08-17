
EVENTS = c("HD","PR","D0","E","1M","2M","3M","6M","taper_start","taper_end","wk2_pt","mo1_pt","mo3_pt","mo6_pt","cGVHD")
EVENTS_sym = c("Healthy Donor",
                "Pre-Conditioning","Day of Infusion","Engraftment",
                "Taper\nStart","Taper\nEnd","2 Week","1 Month","2 Month","3 Month","6 Month","Failure")

COLORs = c("HD"="#379078", "Healthy Donor"="#379078",
            "REF" = "#69818F", "COMP" = "#C94034",
            "aGVHD-"="#69818F", "aGVHD+"="#C94034",
            "GVHD-"="#69818F", "GVHD+"="#377EB8",
            "Success"="#69818F", "Failure"="#377EB8",
            "DFCI"="#f2a73b", "MCC"="#265a92")



GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 0) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



#####--------------------------------------------------------
##### METAGENOMICS
#####--------------------------------------------------------
get_superkingdom_counts <- function(SAMP, SUPKING){
    SAMP_t = gsub(".trimmed","-trimmed",SAMP)
    topA <- read.delim(paste0("/workdir/cfrna/alignment/takara_human_V3/output/HSCT_DFBWH/sample_output/",SAMP_t,"/",SAMP_t,".bracken.extended")) %>% 
            filter(superkingdom == SUPKING) %>%
            mutate(sample_id = SAMP)

    return(topA)
}

get_families_counts <- function(SAMP, FAMs){
    SAMP_t = gsub(".trimmed","-trimmed",SAMP)
    tmp = read.delim(paste0("/workdir/cfrna/alignment/takara_human_V3/output/HSCT_DFBWH/sample_output/",SAMP_t,"/",SAMP_t,".bracken.extended")) %>% 
        group_by(family) %>% summarize(num_reads = sum(new_est_reads)) %>% 
        filter(family %in% FAMs) %>% 
        column_to_rownames("family") %>% t() %>% data.frame() %>% 
        mutate(sample_id = SAMP) %>% remove_rownames()
    INT_COLS = c(FAMs,"sample_id")
    for (i in setdiff(INT_COLS,colnames(tmp))){
        tmp[,i] = 0
    }
    return(tmp)
}

#####--------------------------------------------------------
##### LONGITUIDNAL BOXPLOTS
#####--------------------------------------------------------
make_long_boxplot <- function(meta_plt, MEASUREMENT, y1, y2, y3, YAXIS_LAB=NA, TITLE = NA, EXP_LIM = NA){

    ###-----------------------
    ### get significance
    ###-----------------------

    ### comparisons between groups for each timepoint
    pat_comp = meta_plt %>% 
        filter(!(sample_group %in% c("HD"))) %>%
        group_by(event_sym) %>%
        mutate(MESUR = .data[[MEASUREMENT]]) %>% 
        wilcox_test(MESUR ~ sample_group, paired=FALSE)%>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj") %>% 
        add_xy_position(x = "event_sym",dodge = 0.8)

    ### compare each timepoint to healthy donors
    noCOMP_comp = meta_plt %>% 
        filter(sample_group %in% c("REF","HD")) %>% 
        mutate(MESUR = .data[[MEASUREMENT]]) %>% 
        wilcox_test(MESUR ~ event_sym, paired=FALSE, ref.group = "HD")%>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj") %>% 
        add_xy_position(x = "event_sym") %>% 
        mutate(xmin = xmax) 

    ### compare each timepoint to healthy donors
    COMP_comp = meta_plt %>% 
        filter(sample_group %in% c("COMP","HD")) %>% 
         mutate(MESUR = .data[[MEASUREMENT]]) %>% 
        wilcox_test(MESUR ~ event_sym, paired=FALSE, ref.group = "HD")%>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj") %>% 
        add_xy_position(x = "event_sym") %>% 
        mutate(xmin = xmax) 


    ###-----------------------
    ### make plots
    ###-----------------------

    ### patient plot
    pt_plt <-  meta_plt %>% 
        filter(event_sym != "HD") %>%
        mutate(event_sym = factor(event_sym, levels = c("PR","D0","E","1M","2M","3M","6M"))) %>% 
        ggplot(aes(x=event_sym, y=.data[[MEASUREMENT]], color = sample_group))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(position =  position_jitterdodge(jitter.height=0, 
                                                jitter.width=.25, 
                                                seed=42),
                size = 1, alpha = 0.75, stroke=NA)+
        stat_pvalue_manual(pat_comp,label = "p.adj.signif", tip.length = 0, y.position = y1, group = "sample_group", size = 3, hide.ns= TRUE) +
        stat_pvalue_manual(noCOMP_comp,label = "p.adj.signif", remove.bracket = TRUE, y.position = y2, size = 3, color = "#71b2df", hide.ns= TRUE) +
        stat_pvalue_manual(COMP_comp,label = "p.adj.signif", remove.bracket = TRUE, y.position = y3, size = 3, color = "#dca343", hide.ns= TRUE) 
        
    if(!is.na(EXP_LIM)){
        pt_plt <- pt_plt + expand_limits(y=EXP_LIM)
    }

    ### healthy donor plot
    hd_plt <- meta_plt %>% 
        filter(event_sym == "HD") %>%
        ggplot(aes(x=event_sym, y=.data[[MEASUREMENT]], color = sample_group))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(position =  position_jitterdodge(jitter.height=0, 
                                                jitter.width=.25, 
                                                seed=42),
                size = 1, alpha = 0.75, stroke=NA)

    ### match range
    lim <- range(c(layer_scales(pt_plt)$y$range$range, layer_scales(hd_plt)$y$range$range))

    pt_plt <- pt_plt + ylim(lim)
    hd_plt <- hd_plt + ylim(lim)

    ### add aesthetics
    THEMES <- list(theme_prevail(),
                    theme(
                        axis.title = element_blank(),
                        title = element_blank(),
                        legend.position = "none"),
                    labs(y = "Average correlation to\nhealthy reference"),                    theme( 
                        axis.title.y = element_text(size = 8),
                        strip.background = element_blank(),
                        strip.text.x = element_blank()),
                    scale_color_manual(values = c("REF" = "#71b2df", "COMP" = "#dca343","HD" = "#379078")))

    pt_plt <- pt_plt + THEMES + theme(plot.margin = unit(c(5.5,0,5.5,5.5), unit="pt"))
    hd_plt <- hd_plt + THEMES + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),plot.margin = unit(c(5.5,5.5,5.5,0), unit="pt"))

    if(is.na(YAXIS_LAB)){
        pt_plt <- pt_plt + theme(axis.title.y = element_blank())
        hd_plt <- hd_plt + theme(axis.title.y = element_blank())
    }else(
        pt_plt <- pt_plt + labs(y = YAXIS_LAB)

    )

    ### organize into final plot
    final_plt <- grid.arrange(pt_plt, hd_plt, layout_matrix = rbind(c(1,1,1,1,1,2)))

    if(!is.na(TITLE)){
        final_plt <- grid.arrange(pt_plt, hd_plt, 
                            layout_matrix = rbind(c(1,1,1,1,1,2)),
                            top = textGrob(TITLE, gp=gpar(fontsize=8), vjust = 1.5))
    }
    
    return(final_plt)
}


#####--------------------------------------------------------
##### CORRELATION MEASUREMENTS
#####--------------------------------------------------------

calc_samp_avg <- function(SAMP,
                                matrx,
                                mdf,
                                ref_ids,
                                ref_group){

    if (ref_group =="ES_REF"){
        EVSYM = mdf %>% filter(sample_id == SAMP) %>% pull(event_sym)
        eref_ids <- mdf %>% filter(event_sym == EVSYM) %>% 
                            filter(sample_id %in% ref_ids) %>% 
                            filter(sample_id != SAMP) %>% 
                            pull(sample_id)

        avg_ref_df <- matrx[SAMP,eref_ids,drop=FALSE] %>%
                    mutate(sample_id = row.names(.), avg_ref = rowSums(.) / length(eref_ids)) %>%
                    select(sample_id, avg_ref) %>% remove_rownames()
    }

    if (ref_group == "HD"){
        eref_ids = ref_ids[ref_ids != SAMP]

        avg_ref_df <- matrx[SAMP,eref_ids] %>%
                mutate(sample_id = row.names(.), avg_ref = rowSums(.) / length(eref_ids)) %>%
                select(sample_id, avg_ref) %>% remove_rownames()
}

    return(avg_ref_df)
}

## CORRELATION
get_corr_mtx <- function(cnts,
                            genes_to_use = NA,
                            cor_method = "pearson"
                            ){

    if (is.na(genes_to_use)){
        counts_cor <- cor( cnts, method =  cor_method) %>% data.frame()
        return(counts_cor)
    }else{
        counts_cor <- cor( cnts[genes_to_use,], method =  cor_method) %>% data.frame()
        return(counts_cor)
    }
    
}


get_avg_cor <- function(cnts,
                        ref_ids,
                        meta_data,
                        genes_to_use = NA,
                        ref_group = "HD",
                        cor_method = "pearson"){

    # get sample x sample correlation matrix
    cor_mtrx <- get_corr_mtx(cnts,cor_method = cor_method)


    # calculate average correlation to each reference sample
    all_samps <- lapply(unique(meta_data$sample_id), 
                        calc_samp_avg, 
                        cor_mtrx, 
                        meta_data, 
                        ref_ids,
                        ref_group)
    
    cor_ref_df <- do.call("rbind",all_samps) %>% 
                    data.frame() %>% dplyr::rename(.,avg_cor_ref = avg_ref)

    return(cor_ref_df)    
}
