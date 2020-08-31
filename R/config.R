# set config profile
pair = F
if(pair == T){
  sample_varname <- "Sample}ID"
  time_varname <- "Stage"
  pairID_varname <- "ID"

  time_name   <- c("Before", "After")
  input_micro <- "inst/extdata/abundance.profile"
  input_metadata <- "inst/extdata/phenotype.csv"

  # defalut parameter
  ProjectID   <- "MicroPair"
  time_colour <- c("#FB8072","#BEBADA")
  outputDir <- "./MicroPair"

}else{
  group_name <- c("g0", "g1")
  time_colour <- c("#FB8072","#BEBADA")
  group_varname <- "Group"
}

mytheme <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 10,color = "black", face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, color = "black",face = "bold"),
        axis.text = element_text(size = 9, color = "black"), axis.ticks.length = unit(-0.05, "in"),
        axis.text.y = element_text(margin = unit(c(0.3,0.3, 0.3, 0.3), "cm"), size = 9),
        axis.text.x = element_text(margin = unit(c(0.3,0.3, 0.3, 0.3), "cm")),
        text = element_text(size = 8,color = "black"),
        strip.text = element_text(size = 9, color = "black", face = "bold"),
        panel.grid = element_blank())
