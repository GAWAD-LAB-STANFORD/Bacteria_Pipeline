suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
results_folder <- args[1]
sample_coverage_pattern <- args[2]
setwd(results_folder)

coverage_scatter <- function(data) {
  ggplot(data, aes(x=depth, y=breadth)) + ggtitle("Coverage Scatter") +
    geom_point() + facet_wrap(~ sample) + geom_smooth(method=lm, linetype="dashed",color="darkred") +
    scale_y_continuous(expand = c(0, 0)) + #, limits = c(0, 1)) +#geom_text(label=data$chr) +
    theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
          panel.grid.major = element_line(color = "black"), legend.position = "none",#legend.text=element_text(size=15), 
          panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
          panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
  ggsave("fig_coverage_scatter.pdf", width = 11, height = 8.5)
}
coverage_depth_density <- function(data) {
  ggplot(data, aes(x=depth)) + ggtitle("Depth of Coverage - medians labeled") +#ggtitle(sprintf("Depth of Coverage - Median depth: %s", median(data$depth))) + 
    geom_density(aes(x=depth, group=sample, fill=sample), adjust=2, alpha=0.5) +
    facet_wrap(~ sample) + geom_vline(aes(xintercept=sample_median_depth, color=sample), linetype="dashed") +
    scale_x_continuous(expand = c(0, 0)) + geom_text(aes(x = -Inf, y = -Inf, label = sprintf("%s",sample_median_depth)), hjust = -0.1, vjust = -1) +
    theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
          panel.grid.major = element_line(color = "black"), legend.position = "none",#legend.text=element_text(size=15), 
          panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
          panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
  ggsave("fig_coverage_depth_density.pdf", width = 11, height = 8.5)
}
coverage_breadth_density <- function(data) {
  ggplot(data, aes(x=breadth)) + ggtitle("Breadth of Coverage - medians labeled") +#ggtitle(sprintf("Breadth of Coverage - Median fraction: %s", median(data$breadth))) + 
    geom_density(aes(x=breadth, group=sample, fill=sample), adjust=2, alpha=0.5) +
    facet_wrap(~ sample) + geom_vline(aes(xintercept=sample_median_breadth, color=sample), linetype="dashed") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1.0)) +
    geom_text(aes(x = -Inf, y = -Inf, label = sprintf("%s",sample_median_breadth)), hjust = -0.1, vjust = -1) +
    theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
          panel.grid.major = element_line(color = "black"), legend.position = "none",#legend.text=element_text(size=15), 
          panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
          panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
  ggsave("fig_coverage_breadth_density.pdf", width = 11, height = 8.5)
}
summarize_stats <- function(data) {
  medians <- data %>% group_by(sample) %>% summarize(sample_median_breadth = median(breadth))
  data <- left_join(data, medians, by = "sample")
  medians <- data %>% group_by(sample) %>% summarize(sample_median_depth = median(depth))
  data <- left_join(data, medians, by = "sample")
  write.table(data, "coverage_summarized_stats.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  return(data)
}
# Read in all sample files from a directory into single dataframe
coverage_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("chr", "depth", "breadth", "sample"))
file_names <- dir(results_folder, pattern = sample_coverage_pattern)
if (length(file_names)>=1) {
    for (i in 1:length(file_names)) {
    data <- read.table(file_names[i], sep = "\t", header = TRUE, stringsAsFactors=FALSE)
    label <- rep(substr(gsub("_.*", "", file_names[i]), 1, 2), length(data$chr))
    temp <- data.frame(cbind(data$chr, data$covered_features, data$breadth_coverage_fraction, label))
    temp <- setNames(temp, c("chr", "depth", "breadth", "sample"))
    coverage_df <- rbind(coverage_df, temp)
    rm(data, temp)
    }
    coverage_df$depth <- as.numeric(as.character(coverage_df$depth))
    coverage_df$breadth <- as.numeric(as.character(coverage_df$breadth))
    coverage_df <- summarize_stats(coverage_df)
    coverage_depth_density(coverage_df)
    coverage_breadth_density(coverage_df)
    coverage_scatter(coverage_df)
    cat("Graphed coverage")
}