theme_cfrna_slides <- function(){
  theme_bw() %+replace%
    theme(
          panel.border = element_rect(color='black', fill=NA),
          title = element_text(family="Helvetica", size=15, color='black'),
          axis.title = element_text(family="Helvetica", size=15, color='black'),
          axis.ticks = element_line(color='black'),
          axis.text = element_text(family="Helvetica", size=15, color='black'),
          strip.text = element_text(family="Helvetica", size=15, color='black'),
          strip.background = element_rect(fill="white", color="black"),
          legend.position = "none")
}


theme_cfrna_print <- function(){
  theme_bw() %+replace%
    theme(
          panel.border = element_rect(color='black', fill=NA),
          title = element_text(family="Helvetica", size=8, color='black'),
          axis.title = element_text(family="Helvetica", size=8, color='black'),
          axis.ticks = element_line(color='black'),
          axis.text = element_text(family="Helvetica", size=6, color='black'),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank()
          )

}
