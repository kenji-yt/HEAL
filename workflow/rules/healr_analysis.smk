rule make_healr_polyploid_input_dir:
    input:
        final_sorting_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    output:
        directory("results/healr/input_dir/polyploids/{sample}"),
    log:
        "results/logs/healr/making_{sample}_input_dir.log"
    run:
        with open(log[0], "w") as log_file:
            for prog in PROGENITORS:
                bam_types = [("ref", ""), ("unk", "unknown")] if COUNT_UNKNOWN else [("ref", "")]
                # REF reads in <sample>/<prog>/, unknown reads in <sample>/<prog>/unknown/
                for bam_type, sub_dir in bam_types:
                    bam = f"results/eagle_rc/{wildcards.sample}/{wildcards.sample}_classified_{prog}.{bam_type}.bam"
                    if not os.path.isfile(bam):
                        raise ValueError(f"Expected EAGLE-RC output {bam} not found.")
                    link_dir = os.path.join(output[0], prog, sub_dir)
                    os.makedirs(link_dir, exist_ok=True)
                    link = os.path.join(link_dir, os.path.basename(bam))
                    os.symlink(os.path.relpath(bam, link_dir), link)  # relative link, valid if results/ is moved as a whole
                    log_file.write(f"{link} -> {bam}\n")

def get_all_healr_inputs(wildcards):

    sample_path = "results/healr/input_dir/polyploids/"
    sample_dict = {sample: sample_path + sample for sample in SAMPLES}

    progenitor_path = "results/healr/input_dir/subgenomes/"
    gc_dict = {progenitor + "_gc": progenitor_path + progenitor + "/" + progenitor + f"_{BIN_SIZE}_gc.bed" for progenitor in PROGENITORS}
    
    mappability_path = "results/genmap/"
    map_dict = {progenitor + "_map": mappability_path + progenitor + "/index_dir" for progenitor in PROGENITORS}

    merged_dict = sample_dict | gc_dict | map_dict
    return list(merged_dict.values())


rule healr_analysis:
    input:
        get_all_healr_inputs,
        genespace_log="results/logs/genespace/genespace_run.log",
    output:
        stats=directory("results/healr/stats"),
        healr_list=directory("results/healr/healr_list"),
        figures=directory("results/healr/figures"),
        filtering_log="results/healr/filtering_log.txt",
    log:
        "results/logs/healr/analysis.log"
    threads:workflow.cores
    conda:
        "../envs/healr.yaml"
    params:
        input_dir="results/healr/input_dir",
        is_paired=f"{IS_PAIRED}",
        genespace_dir="results/genespace/run_dir",
        data_type=DATA_TYPE,
        healr_version=HEALR_VERSION,
        count_unknown=COUNT_UNKNOWN
    script:
        f"{workflow.basedir}/scripts/healr_analysis.R" 