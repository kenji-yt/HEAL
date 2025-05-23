sink(snakemake@log[[1]], split=TRUE, append=TRUE)

library(devtools)

devtools::install_github("kenji-yt/healr", quiet = TRUE)

library(healr)

input_directory <- snakemake@params[["input_dir"]]
is_paired <- as.logical(snakemake@params[["is_paired"]])
genespace_directory <- snakemake@params[["genespace_dir"]]
num_threads <- snakemake@threads

plot_out_dir <- paste0(dirname(input_directory),"/figures")
stats_out_dir <- paste0(dirname(input_directory),"/stats")


count_list <- count_heal_data(input_dir = input_directory,
                n_threads= num_threads,
                paired_end = is_paired,
                full_output = FALSE)

# Add options in config here
filt_list <- filter_bins(count_list)

# If DNA then correct with GC: or not

cn_list <- get_copy_number(heal_list = filt_list,
                n_threads = num_threads)

alignment_list <- get_heal_alignment(heal_list = cn_list,
                                     genespace_dir = genespace_directory,
                                     n_threads = num_threads)

plot_alignment(heal_list = cn_list,
               alignment = alignment_list,
               output_dir = plot_out_dir,
               n_threads = num_threads)

plot_heal_heat_map(alignment = alignment_list,
                   output_dir = plot_out_dir)

summary_aln <- summarize_aln(alignment = alignment_list,
                             n_threads = num_threads)

write_aln_summary(alignment_summary = summary_aln,
                 output_dir = stats_out_dir)

