rule make_healr_input_dir:
    input:
        snake-genespace-rep-report=="results/snake-genespace-reproducibility_report.txt"
        snake-EAGLE-RC-rep-report="results/snake-EAGLE-RC-reproducibility_report.txt"
        genmap_log="results/logs/bin_information/genmap_{progenitor}.log",
        gc_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_gc.bed",
        bin_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
        healr_download_log="results/logs/healr_download/download.log"
    output:
        directory("results/healr_input_dir")
