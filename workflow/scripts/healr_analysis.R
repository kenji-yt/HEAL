# Open file connection for writing
log <- file(snakemake@log[[1]], open = "wt") 

sink(log, split=TRUE, append=TRUE)
sink(log, append=TRUE, type="message")

# Install and load healr
healr_ref <- paste0("v", snakemake@params[["healr_version"]])
if (!requireNamespace("healr", quietly = TRUE) ||
    !identical(utils::packageDescription("healr")$RemoteRef, healr_ref))
    remotes::install_github(paste0("kenji-yt/healr@", healr_ref), quiet = TRUE, upgrade = "never", force = TRUE)

library(healr)

# Get parameters
input_directory <- snakemake@params[["input_dir"]]
is_paired <- as.logical(snakemake@params[["is_paired"]])
genespace_directory <- snakemake@params[["genespace_dir"]]
data_type <- snakemake@params[["data_type"]]
count_unknown <- as.logical(snakemake@params[["count_unknown"]])
num_threads <- snakemake@threads

plot_out_dir <- paste0(dirname(input_directory),"/figures")
stats_out_dir <- paste0(dirname(input_directory),"/stats")

# Analyses
count_list <- count_heal_data(input_dir = input_directory,
                n_threads= num_threads,
                paired_end = is_paired,
                full_output = FALSE,
                unknown_reads = count_unknown)


healr_lists_out_dir <- paste0(dirname(input_directory),"/healr_list")

write_heal_list(heal_list = count_list,
                output_dir = healr_lists_out_dir)


filt_list <- filter_bins(count_list, log_file = paste0(dirname(input_directory),"/filtering_log.txt"))

if(data_type == "DNA"){
    dir.create(plot_out_dir)
    filt_list <- correct_gc(heal_list = filt_list,
               n_threads = num_threads, 
               output_dir = plot_out_dir)
}

cn_list <- get_cn(heal_list = filt_list,
                n_threads = num_threads)

alignment_list <- get_heal_aln(heal_list = cn_list,
                                     genespace_dir = genespace_directory,
                                     n_threads = num_threads)

plot_aln(heal_list = cn_list,
               alignment = alignment_list,
               output_dir = plot_out_dir,
               n_threads = num_threads,
               add_bins = "all")

plot_all_bins(heal_list = cn_list,
               plot_cn = TRUE,
               output_dir = plot_out_dir)

plot_riparian(heal_list = cn_list,
               alignment = alignment_list,
               genespace_dir = genespace_directory,
               output_dir = plot_out_dir)

plot_aln_heat(alignment = alignment_list,
                   output_dir = plot_out_dir)

summary_aln <- summarize_aln(alignment = alignment_list,
                             n_threads = num_threads)

write_aln_summary(alignment_summary = summary_aln,
                 output_dir = stats_out_dir)