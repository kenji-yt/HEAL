sink(snakemake@log[[1]], split=TRUE, append=TRUE)
sink(snakemake@log[[1]], split=TRUE, append=TRUE, type="message")

library(devtools)

devtools::install_github("kenji-yt/healr", quiet = TRUE)

library(healr)

input_directory <- snakemake@params[["input_dir"]]
is_paired <- as.logical(snakemake@params[["is_paired"]])
genespace_directory <- snakemake@params[["genespace_dir"]]
data_type <- snakemake@params[["data_type"]]
num_threads <- snakemake@threads

plot_out_dir <- paste0(dirname(input_directory),"/figures")
stats_out_dir <- paste0(dirname(input_directory),"/stats")

    
count_list <- count_heal_data(input_dir = input_directory,
                n_threads= num_threads,
                paired_end = is_paired,
                full_output = FALSE)


healr_lists_out_dir <- paste0(dirname(input_directory),"/healr_list")

write_heal_list(heal_list = count_list,
                output_dir = healr_lists_out_dir)


filt_list <- filter_bins(count_list, log_file = paste0(dirname(input_directory),"/filtering_log.txt"))

if(data_type == "DNA"){
    filt_list <- correct_gc(heal_list = filt_list,
               n_threads = num_threads, 
               output_dir = plot_out_dir)
}

cn_list <- get_copy_number(heal_list = filt_list,
                n_threads = num_threads)

alignment_list <- get_heal_alignment(heal_list = cn_list,
                                     genespace_dir = genespace_directory,
                                     n_threads = num_threads)

plot_alignment(heal_list = cn_list,
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

plot_heal_heat_map(alignment = alignment_list,
                   output_dir = plot_out_dir)

summary_aln <- summarize_aln(alignment = alignment_list,
                             n_threads = num_threads)

write_aln_summary(alignment_summary = summary_aln,
                 output_dir = stats_out_dir)

