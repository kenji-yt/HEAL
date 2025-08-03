rule make_healr_polyploid_input_dir:
    input:
        final_sorting_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    output:
        directory("results/healr/input_dir/polyploids/{sample}"),
    log:
        "results/logs/healr/making_{sample}_input_dir.log"
    params:
        progenitor_list=" ".join(PROGENITORS)
    shell:
        '''
        for prog in {params.progenitor_list}; do
            out_dir="results/healr/input_dir/polyploids/{wildcards.sample}/${{prog}}/"
            mkdir -p "${{out_dir}}" 2> "{log}"
            ref_bam=$(find results/eagle_rc/{wildcards.sample}/ -name *_classified_${{prog}}.ref.bam)
            ln -s "../../../../../../${{ref_bam}}" "${{out_dir}}" 2> "{log}"
        done
        '''

def get_all_healr_inputs(wildcards):

    sample_path = "results/healr/input_dir/polyploids/"
    sample_dict = {sample: sample_path + sample for sample in SAMPLES}

    progenitor_path = "results/healr/input_dir/progenitors/"
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
        directory("results/healr/stats")
    log:
        "results/logs/healr/analysis.log"
    threads:workflow.cores
    params:
        input_dir="results/healr/input_dir",
        is_paired=f"{IS_PAIRED}",
        genespace_dir="results/genespace/run_dir",
    script:
        f"{workflow.basedir}/scripts/healr_analysis.R" 