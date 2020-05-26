suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(reshape2)
  library(directlabels)
  library(rjson)
  library(treeio)
  library(ggtree)
})

args <- commandArgs(trailingOnly = TRUE)
project <- args[1]
sample_read_count_filename <- args[2]
summed_read_targets_filename <- args[3]
contig_read_targets_filename <- args[3]
contig_data_filename <- args[4]
kraken_db_types <- unlist(strsplit(args[5], "-"))
kraken_jtree_suffix <- args[6]
blast_db_types <- unlist(strsplit(args[7], "-"))
blast_results_suffix <- args[8]
ncbi_annotations_dir <- args[9]


# Common themes and functions ------------------------------------------------------
ggplot_theme <- theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
                      panel.grid.major = element_line(color = "black"), legend.text=element_text(size=15), 
                      panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
                      panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
ggplot_theme_small_legend <- theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
                                   panel.grid.major = element_line(color = "black"), legend.title = element_text(size = 5), legend.text = element_text(size = 5), 
                                   panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
                                   panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
ggplot_theme_no_legend <- theme(axis.line.y = element_line(size=.1,color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
                                panel.grid.major = element_line(color = "black"), legend.position = "none",
                                panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
                                panel.grid.major.y = element_line( size=.1, color="black" ), plot.title = element_text(size=20))
ggtree_theme <- theme(axis.line.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=10, lineheight=0.2, color="black"),
                      panel.grid.major = element_line(color = "black"), legend.text=element_text(size=15), 
                      panel.background = element_rect(fill="white"), panel.grid.major.x = element_blank() , 
                      panel.grid.major.y = element_blank(), plot.title = element_text(size=20))
grid_table_theme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.9)), colhead = list(fg_params=list(cex = 1.0)))

parse_long_sample_name <- function(df, filename = NULL) {
  if (str_count(df$sample, "_")[1] == 6) {
    df <- df %>%
      separate(sample, sep = "_", c("project", "cell_type", "tissue_origin", "sequencing_type", "TBID", "sample", "sample_number")) %>%
      mutate(sample_number = str_remove(sample_number, "S")) 
    if (!is.null(filename)) {
      write.table(df, filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  }
  return(df)
}


# Sample read counts ------------------------------------------------------
sample_read_counts_df <- read.table(sample_read_count_filename, sep="\t", stringsAsFactors = FALSE, header = TRUE)
sample_read_counts_df <- parse_long_sample_name(sample_read_counts_df, sample_read_count_filename)
plot1 <- ggplot(sample_read_counts_df, aes(x = 1, y = read_count)) + 
  geom_boxplot() + geom_jitter(aes(color = sample), height = 0) + 
  labs(x = "", y = "Read count", title  = "Sample read counts") +
  ggplot_theme + xlim(c(0,2))
pdf(sprintf("%s.sample_read_counts.pdf", project), width = 15, height = 8.5)
grid.arrange(tableGrob(sample_read_counts_df[,c("sample", "read_count")], rows = NULL, theme = grid_table_theme), plot1, ncol = 2)
dev.off()

expected_coverage_df <- data.frame(genome_size = 0:10) %>% 
  mutate(mean_cov = (median(sample_read_counts_df$read_count) * 140) / (genome_size * 1e6))
plot1 <- ggplot(expected_coverage_df, aes(x = genome_size, y = mean_cov)) + geom_line() + ggplot_theme +
  labs(x = "Genome size (MB)", y = "Mean coverage depth", title  = "Expected mean coverage depth given sample genome sizes")
ggsave(sprintf("%s.fig_expected_mean_coverage.pdf", project), plot = plot1, width = 11, height = 8.5)
cat("Plotted sample read counts for all samples\n")
rm(sample_read_counts_df, expected_coverage_df, plot1)


# Summed read targets ------------------------------------------------------------
summed_read_targets_df <- read.table(summed_read_targets_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
summed_read_targets_df <- parse_long_sample_name(summed_read_targets_df, summed_read_targets_filename)
summed_read_targets_df <- summed_read_targets_df %>%
  select(sample, human_aligned, contig_aligned, unaligned) %>%
  melt(id.vars = c("sample"), value.name = "read_count", variable.name = "alignment") %>%
  group_by(sample) %>%
  mutate(percent_read_count = read_count / sum(read_count) * 100) %>%
  ungroup()

# Read counts for three different categories: contig, human, and undetermined
#   Read counts from BWA ALN aligning of all reads against their Contigs and Hg38 reference
plot1 <- ggplot(summed_read_targets_df, aes(sample, read_count, fill = reorder(alignment, read_count))) + 
  geom_bar(stat = "identity", position = position_dodge()) + geom_text(aes(label = round(percent_read_count, 2)), vjust=-1, position = position_dodge(0.9)) +
  labs(title = "Human contamination by read count", x = "Sample", y = "Read count", fill = "Read categories") + ggplot_theme
ggsave(sprintf("%s.fig_human_contamination.pdf", project), plot = plot1, width = 11, height = 8.5)
cat("Plotted human contamination for all samples\n")
rm(summed_read_targets_df, plot1)


# Contig read targets ------------------------------------------------------------
contig_read_targets_df <- read.table(contig_read_targets_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (! "read_count" %in% colnames(contig_read_targets_df)) { 
  contig_read_targets_df <- contig_read_targets_df %>%
    group_by(sample, target) %>% 
    summarize(read_count = n()) %>%
    ungroup()
  contig_read_targets_df <- parse_long_sample_name(contig_read_targets_df, contig_read_targets_filename)
  cat("Summarized read targets for all samples\n")
}


# Contig data -------------------------------------------------------------
contig_data_df <- read.table(contig_data_filename, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
if ("fasta_header" %in% colnames(contig_data_df)) {
  contig_data_df <- parse_long_sample_name(contig_data_df)
  contig_data_df <- contig_data_df %>%
  mutate(fasta_header = str_remove(fasta_header, ">")) %>%
  rename(contig = fasta_header) %>%
  mutate(temp_list = sapply(strsplit(contig, "_"), function(.split_vec) {return(list(as.numeric(c(.split_vec[2], .split_vec[4], .split_vec[6]))))}),
              contig_rank = sapply(temp_list, function(.x) {.x[1]}),
              contig_length = sapply(temp_list, function(.x) {.x[2]}),
              contig_graph_cov = sapply(temp_list, function(.x) {.x[3]})
  ) %>% select(-temp_list)
  contig_data_df <- left_join(contig_data_df, contig_read_targets_df[, c("sample", "target", "read_count")], by = c("sample", "contig" = "target"))
  rm(contig_read_targets_df)
  if ("project" %in% colnames(contig_data_df)) {
    contig_data_df <- select(contig_data_df, -project, -cell_type, -tissue_origin, -sequencing_type, -TBID, -sample_number)
  }
  write.table(contig_data_df, contig_data_filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("Summarized contig data for all samples\n")
}

# Contig length, rank, and SPAdes graph coverage across all Samples
contig_data_df <- contig_data_df %>% group_by(sample) %>% mutate(cumsum_contig_length = cumsum(contig_length)) %>% ungroup()
cum_sum_contig_length_plot <- ggplot(mutate(contig_data_df, trunc_contig_graph_cov = ifelse(contig_graph_cov > 100, 100, contig_graph_cov)), 
                                     aes(x = contig_rank, y = cumsum_contig_length, size=(trunc_contig_graph_cov), color = sample)) + 
  geom_line() + facet_grid(. ~ sample) + labs(x = "Contig rank", y = "Contig length") + ggplot_theme_no_legend
cum_sum_contig_cov_plot <- ggplot(mutate(contig_data_df, trunc_contig_graph_cov = ifelse(contig_graph_cov > 100, 100, contig_graph_cov)), 
                                  aes(x = contig_rank, y = trunc_contig_graph_cov, color = sample)) + 
  geom_point() + facet_grid(. ~ sample) + labs(x = "Contig rank", y = "SPAdes graph coverage") + ggplot_theme_no_legend
pdf(sprintf("%s.fig_contig_rank_length_coverage.pdf", project), width = 11, height = 8.5)
grid.arrange(cum_sum_contig_length_plot, cum_sum_contig_cov_plot, nrow = 2)
dev.off()
cat("Plotted contig data for all samples\n")
contig_data_df <- select(contig_data_df, -cumsum_contig_length)
rm(cum_sum_contig_length_plot, cum_sum_contig_cov_plot)


# Functions for BLAST and Kraken2 results ---------------------------------------------
preprocess_blast_results <- function(sample_df, contig_data_df) {
  sample_df <- parse_long_sample_name(sample_df)
  if ("species" %in% colnames(sample_df)) {
    sample_df <- sample_df %>%
      filter(hit_rank == 1) %>% 
      mutate(hit_taxid = factor(hit_taxid), genus = sub(pattern = " .+", "", species))
  }
  sample_df <- filter(sample_df, hit_rank == 1) %>% 
    mutate(tophit_aln_query_fraction = top_hsp_align_len / contig_length) %>%
    group_by(sample) %>% 
    mutate(query_length_rank = 1:length(contig_length)) %>% 
    ungroup()
  sample_df <- left_join(sample_df, contig_data_df, by = c("sample", "contig", "contig_length"))
  return(sample_df)
}
contig_to_blast_alignment_plot <- function(sample_df, blast_type_plot_string) {
  blast_type_plot_string <- tolower(blast_type_plot_string)
  plot1 <- ggplot(data = sample_df, aes(x = tophit_aln_query_fraction, fill = sample)) + 
    facet_wrap(~ sample) + geom_histogram(binwidth = 0.2) +
    labs(title = sprintf("Fraction of contig aligning to top %s BLAST hit", blast_type_plot_string), 
         x = "Fraction of contig aligning to top BLAST hit", y = "Number of contigs") +
    ggplot_theme_no_legend
  return(plot1)
}
heatmap_for_all_samples <- function(blast_df, taxonomic_plot_and_variable_string) {
  taxonomic_variable_string <- tolower(taxonomic_plot_and_variable_string)
  unique_taxonomy <- sort(unique(pull(blast_df, taxonomic_variable_string)))
  sample_to_all_counts_list <- lapply(unique(blast_df$sample), function(.sample) {
    sample_df <- filter(blast_df, sample == .sample) %>% 
      group_by(.dots = taxonomic_variable_string) %>% 
      summarize(total_contig_length = sum(contig_length))
    temp_df <- data.frame(sample = .sample, hit_taxon = unique_taxonomy, stringsAsFactors = FALSE)
    temp_with_counts_df <- left_join(temp_df, sample_df, by = c("hit_taxon" = taxonomic_variable_string)) %>% 
      mutate(total_contig_length = ifelse(is.na(total_contig_length), 0, total_contig_length))
  })
  
  sample_to_all_counts_df <- do.call(rbind, sample_to_all_counts_list)
  plot1 <- ggplot(sample_to_all_counts_df, aes(y = sample, x = hit_taxon)) + 
    geom_tile(aes(fill = total_contig_length)) + ggplot_theme_no_legend +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
    labs(fill ="Total contig length", x = taxonomic_plot_and_variable_string, y = "Sample", title = sprintf("Sum of contig lengths matching %s", taxonomic_variable_string)) + 
    theme(plot.title = element_text(hjust = 0.5))
  return(plot1)
}
taxonomic_color_vector <- function(sample_df, taxonomic_variable_string) {
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  unique_taxonomy <- sort(unique(pull(sample_df, taxonomic_variable_string)))
  taxonomy_colors <- gg_color_hue(length(unique_taxonomy))
  names(taxonomy_colors) <- unique_taxonomy
  return(taxonomy_colors)
}
summarize_sample_metrics <- function(sample_df, taxonomic_variable_string) {
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  summed_sample_df <- sample_df %>%
    mutate(reverse_ranking = (max(contig_rank) + 1) - contig_rank) %>%
    mutate(total_read_count = sum(read_count),
           total_contigs = length(contig), 
           total_contig_length = sum(contig_length),
           total_blast_length = sum(top_hsp_align_len),
           total_contig_graph_cov = sum(contig_graph_cov),
           total_reverse_ranking = sum(reverse_ranking),
           total_blast_alignment_fractions = sum(tophit_aln_query_fraction)) %>%
    group_by(.dots = taxonomic_variable_string) %>%
    mutate(max_reference_len = max(reference_len),
           summed_read_count = sum(read_count),
           summed_contigs = n(),
           summed_contig_length = sum(contig_length),
           summed_blast_length = sum(top_hsp_align_len),
           summed_contig_graph_cov = sum(contig_graph_cov),
           summed_reverse_rankings = sum(reverse_ranking),
           summed_blast_alignment_fractions = sum(tophit_aln_query_fraction)) %>%
    mutate(percent_read_count = summed_read_count / total_read_count * 100,
           percent_contigs = summed_contigs / total_contigs * 100,
           percent_contig_length = summed_contig_length / total_contig_length * 100,
           percent_blast_length = summed_blast_length / total_blast_length * 100,
           percent_contig_graph_cov = summed_contig_graph_cov / total_contig_graph_cov * 100,
           percent_reverse_ranking = summed_reverse_rankings / total_reverse_ranking * 100,
           percent_blast_alignment_fraction = summed_blast_alignment_fractions / total_blast_alignment_fractions * 100,
           percent_reference_covered_by_contigs = summed_contig_length / max_reference_len * 100,
           percent_reference_covered_by_blasts = summed_blast_length / max_reference_len * 100) %>%
    select(sample, taxonomic_variable_string, summed_read_count, summed_contigs, summed_contig_length, summed_blast_length, summed_contig_graph_cov, summed_blast_alignment_fractions,
           percent_read_count, percent_contigs, percent_contig_length, percent_contig_graph_cov, percent_reverse_ranking, 
           percent_blast_alignment_fraction, percent_reference_covered_by_contigs, percent_reference_covered_by_blasts, max_reference_len) %>%
    distinct() %>%
    ungroup()
  return(summed_sample_df)
}
kraken_ggtree_plot <- function(sample_string, db_type) {
  json <- fromJSON(file = sprintf("%s.%s.%s%s", project, sample_string, db_type, kraken_jtree_suffix))
  width <- json$metadata$max_depth + 1
  tree <- read.jtree(sprintf("%s.%s.%s%s", project, sample_string, db_type, kraken_jtree_suffix))
  tree1 <- ggtree(tree, branch.length='none', aes(color=percent_fragments_covered), size = 1) + 
    geom_label(aes(x=branch, label=label), vjust=-1) + geom_label(aes(x=branch, label=percent_fragments_covered)) +
    geom_tiplab(size=5, color="black") + 
    labs(title = "Percent coverage from kraken2 results", x = "Depth of identification") +
    ylim(0, width/2) + xlim(0, width) + theme(legend.position="bottom") + ggtree_theme +
    scale_color_continuous(low='red', high='royalblue1')
  return(tree1)
}
metric_bar_plots_list <- function(summed_sample_df, taxonomic_plot_string, taxonomic_variable_string, taxonomic_color_vector) {
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  summed_sample_df <- summed_sample_df %>% top_n(10, summed_contig_length)
  plot1 <- ggplot(summed_sample_df, aes(reorder(get(taxonomic_variable_string), percent_read_count), percent_read_count, fill = get(taxonomic_variable_string))) +
    geom_bar(stat="identity") + geom_text(aes(label=round(summed_read_count,2)), vjust=-1) +
    labs( x = taxonomic_plot_string, y = "Percent of all reads", title = sprintf("%s read count", taxonomic_plot_string)) + 
    ylim(0, 100) + ggplot_theme_no_legend + scale_color_manual(values = taxonomic_color_vector)
  plot2 <- ggplot(summed_sample_df, aes(reorder(get(taxonomic_variable_string), percent_contig_length), percent_contig_length, fill = get(taxonomic_variable_string))) +
    geom_bar(stat="identity") + geom_text(aes(label=round(summed_contig_length,2)), vjust=-1) +
    labs( x = taxonomic_plot_string, y = "Percent of all contig lengths", title = sprintf("%s contig length", taxonomic_plot_string)) + 
    ylim(0, 100) + ggplot_theme_no_legend + scale_color_manual(values = taxonomic_color_vector)
  plot3 <- ggplot(summed_sample_df, aes(reorder(get(taxonomic_variable_string), percent_reference_covered_by_blasts), percent_reference_covered_by_blasts, fill = get(taxonomic_variable_string))) +
    geom_bar(stat="identity") + geom_text(aes(label=round(percent_reference_covered_by_blasts,2)), vjust=-1) +
    labs( x = taxonomic_plot_string, y = "Percent of reference genome length", title = sprintf("%s contig aligned length\n compared to reference length", taxonomic_plot_string)) + 
    ylim(0, 100) + ggplot_theme_no_legend + scale_color_manual(values = taxonomic_color_vector)
  return_list = list(plot1, plot2, plot3)
  return(return_list)
}
taxonomic_proportion_of_contig_length_plot <- function(sample_df, taxonomic_plot_string, taxonomic_variable_string, taxonomic_color_vector) {
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  cumprop_df <- do.call(rbind, lapply(1:nrow(sample_df), function(.row_num) {
    data.frame(table(pull(sample_df, taxonomic_variable_string)[1:.row_num])) %>% 
      mutate(prop_taxonomy = Freq / sum(Freq), contig_length = sum(sample_df$contig_length[1:.row_num]))
    })
  )
  plot1 <- ggplot(cumprop_df, aes(contig_length, y = prop_taxonomy, color = Var1)) + 
    geom_line() + ggplot_theme + scale_color_manual(values = taxonomic_color_vector) + ggplot_theme_no_legend +
    labs(x = "Cumulative contig length", y = "Proportion of cumulative length", title = sprintf("%s contig length proportions", taxonomic_plot_string)) +
    geom_dl(aes(label = Var1),  method = list("last.points", cex = 0.8, hjust=1.2, vjust=1.2))
  return(plot1)
}
scaled_taxonomic_proportion_of_contig_length_plot <- function(sample_df, taxonomic_plot_string, taxonomic_variable_string, taxonomic_color_vector) {
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  cumprop_df <- do.call(rbind, lapply(1:nrow(sample_df), function(.row_num) {
    sample_df[1:.row_num,] %>% group_by(.dots = taxonomic_variable_string) %>% 
      summarize(cum_taxonomy_contig_length = sum(contig_length)) %>% 
      mutate(total_cumsum_contig_length = sum(cum_taxonomy_contig_length), prop_cum_taxonomy_contig_length = cum_taxonomy_contig_length / total_cumsum_contig_length)
    })
  )
  plot1 <- ggplot(cumprop_df, aes(x = total_cumsum_contig_length, y = prop_cum_taxonomy_contig_length , color = get(taxonomic_variable_string))) + 
    geom_line() + ggplot_theme + scale_color_manual(values = taxonomic_color_vector) + ggplot_theme_no_legend +
    labs(x = "Cumulative contig length", y = "Proportion of cumulative length", title = sprintf("%s contig length proportions, scaled by group", taxonomic_plot_string)) +
    geom_dl(aes(label = get(taxonomic_variable_string)), method = list("last.points", cex = 1, hjust=1.2, vjust=1.2))
  return(plot1)
}
export_with_ncbi_annotations <- function(summed_df, ncbi_annotation_string, taxonomic_variable_string, blast_db_string) {
  ncbi_annotation_string <- tolower(ncbi_annotation_string)
  taxonomic_variable_string <- tolower(taxonomic_variable_string)
  if (!is.na(ncbi_annotations_dir)) {
    annotation_df <- read.table(sprintf("%s/%s_annotations.tsv", ncbi_annotations_dir, ncbi_annotation_string), sep="\t", header = TRUE, quote="")
    summed_df <- plyr::join(summed_df, annotation_df, by = c(taxonomic_variable_string), type = "left", match = "first")
  }
  if (ncbi_annotation_string == "bacteria") {
    write.table(summed_df, sprintf("%s.summarized_%s_%s.tsv", project, blast_db_string, taxonomic_variable_string), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    write.table(summed_df, sprintf("%s.summarized_%s.tsv", project, blast_db_string), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  cat(sprintf("Summarized and annotated %s BLAST results for all samples\n", blast_db_string))
}


# Nucleotide, plasmid, and viral BLAST results ------------------------------------------------
taxonomic_variable_string_list <- c("Genus", "Plasmid", "Phage")
ncbi_annotation_string_list <- c("Bacteria", "Plasmid", "Viral")

# Load and preprocess all BLAST result dataframes
# Plot contig alignment bar plots, taxonomy heatmaps, 
# and get list of colors corresponding to taxonomy variables for plot colors later
blast_results_df_list <- list()
alignments <- list()
heatmaps <- list()
taxonomic_color_list <- list()
for (i in 1:length(blast_db_types)) {
  filename <- sprintf("%s.%s%s", project, blast_db_types[i], blast_results_suffix)
  if (file.exists(filename)) {
    df <- read.table(filename, sep="\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
    all_alignments_df <- preprocess_blast_results(df, contig_data_df)
    df <- filter(all_alignments_df, tophit_aln_query_fraction >= 0.9)
    if (nrow(df) > 0) {
      blast_results_df_list <- c(blast_results_df_list, list(df))
      heatmaps <- c(heatmaps, list(heatmap_for_all_samples(blast_results_df_list[[i]], taxonomic_variable_string_list[i])))
      taxonomic_color_list <- c(taxonomic_color_list, list(taxonomic_color_vector(blast_results_df_list[[i]], taxonomic_variable_string_list[i])))
    }
    if (nrow(all_alignments_df) > 0) {
      alignments <- c(alignments, list(contig_to_blast_alignment_plot(all_alignments_df, blast_db_types[i])))
    }
  } else {
    cat(sprintf("Missing %s\n", filename))
  }
}
unique_samples <- unique(contig_data_df$sample)
rm(contig_data_df, df, all_alignments_df)

# End program if no BLAST results found
if (length(blast_results_df_list) < 1) {
  cat("No BLAST results found, ending script\n")
  stop()
}

pdf(sprintf("%s.fig_contig_alignments_to_blast_hits.pdf", project), width = 8.5, height = 11)
grid.arrange(grobs = alignments, ncol = 1)
dev.off()
cat("Plotted contig alignments to BLAST hits for all samples\n")

pdf(sprintf("%s.fig_blast_results_heatmaps.pdf", project), width = 8.5, height = 11)
grid.arrange(grobs = heatmaps, ncol = 1)
dev.off()
cat("Plotted heatmaps of BLAST results for all samples\n")
rm(alignments, heatmaps)


# Individual grid figures for each Sample with 
#   kraken phylogenetic tree plot
#   plots for summarized metrics of nucleotide Genus, plasmid Plasmids, and viral Phages
# Consolidate summarized metric data into a list of dataframes
summarized_df_list <- list(data.frame(), data.frame(), data.frame())
summarized_nucleotide_species_df <- data.frame()
cat(sprintf("Plotting kraken and BLAST results for %s samples\n", length(unique_samples)))
count <- 1
for (current_sample in unique_samples) {
  grid_list <- list()
  grid_layout <- data.frame()
  variable_height <- 0
  
  for (i in 1:length(blast_results_df_list)) {
    sample_df <- blast_results_df_list[[i]] %>%
      filter(sample == current_sample)
    if (nrow(sample_df) > 0) {
      summed_sample_df <- summarize_sample_metrics(sample_df, taxonomic_variable_string_list[i])
      
      grid_list <- c(grid_list, list(kraken_ggtree_plot(current_sample, kraken_db_types[i])))
      grid_list <- c(grid_list, metric_bar_plots_list(summed_sample_df, ncbi_annotation_string_list[i], taxonomic_variable_string_list[i], taxonomic_color_list[[i]]))
      grid_list <- c(grid_list, list(scaled_taxonomic_proportion_of_contig_length_plot(sample_df, ncbi_annotation_string_list[i], taxonomic_variable_string_list[i], taxonomic_color_list[[i]])))
      
      summarized_df_list[[i]] <- rbind(summarized_df_list[[i]], summed_sample_df)
      if (taxonomic_variable_string_list[i] == "Genus") {
        summed_sample_df <- summarize_sample_metrics(sample_df, "Species")
        summarized_nucleotide_species_df <- rbind(summarized_nucleotide_species_df, summed_sample_df)
      }
      
      variable_height = variable_height + 8
      lay_start <- (nrow(grid_layout)*5)+1
      new_lay <- c(lay_start, lay_start, lay_start+1, lay_start+2, lay_start+3, lay_start+4, lay_start+4)
      if (nrow(grid_layout) == 0) {
        grid_layout <- rbind(new_lay)
      } else {
        grid_layout <- rbind(grid_layout, new_lay)
      }
    }
  }
  
  pdf(sprintf("%s.%s.fig_kraken_blast_results.pdf", project, current_sample), width = 30, height = variable_height)
  grid.arrange(grobs = grid_list, layout_matrix = grid_layout)
  dev.off()
  cat(sprintf("\t%s - %s\n", count, current_sample))
  count = count + 1
}
cat("\tDone\n")

# Annotate all consolidated summarized taxonomy metrics and export as data table files
for (i in 1:length(summarized_df_list)) {
  if (nrow(summarized_df_list[[i]]) > 0) {
    export_with_ncbi_annotations(summarized_df_list[[i]], ncbi_annotation_string_list[i], taxonomic_variable_string_list[i], blast_db_types[i])
    if (taxonomic_variable_string_list[i] == "Genus") {
      export_with_ncbi_annotations(summarized_nucleotide_species_df, ncbi_annotation_string_list[i], "Species", blast_db_types[i])
    }
  }
}
