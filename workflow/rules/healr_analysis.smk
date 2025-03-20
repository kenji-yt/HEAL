rule make_healr_polyploid_input_dir:
    input:
        final_sorting_log="results/logs/eagle_rc/restoring_chr_names/{sample}.log",
    output:
        directory("results/healr/input_dir/polyploids/{sample}")
    log:
        "results/logs/healr/making_input_dir.log"
    params:
        progenitor_list="{PROGENITORS}"
    shell:
        '''
        for prog in {params.progenitor_list}, do
            out_dir='results/healr/input_dir/polyploids/{wildcards.sample}/${prog}/'
            mkdir -p ${out_dir}
            ref_bam=$(find results/eagle-rc/{wildcards.sample}/ -name *_classified_${prog}.ref.bam)
            ln -s ${ref_bam} ${out_dir}
        done
        '''


rule healr_analysis:
    input:
        input_dir_log="results/logs/healr/making_input_dir.log",
        genespace_log="results/logs/genespace/genespace_run.log"
    output:
        directory("results/healr/stats")
    log:
        "results/logs/healr/analysis.log"
    threads:workflow.cores
    params:
        input_dir="results/healr/input_dir",
        is_paired=f"{IS_PAIRED}"
        genespace_dir="results/genespace/run_dir"
    script:
        "scripts/healr_analysis.R" 