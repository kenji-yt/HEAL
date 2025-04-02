rule define_bins:
    input:
        lambda wildcards: EAGLE_RC.get_assembly(wildcards.progenitor),
    output:
        bin_bed="results/healr/input_dir/progenitors/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
    log:
        "results/logs/bin_information/define_bins/{progenitor}.log",
    params:
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell:
        """
        samtools faidx {input} 
        bedtools makewindows -g {input}.fai -w {params.bin_size} > {output.bin_bed}
        """

rule count_gc:
    input:
        assembly=lambda wildcards: EAGLE_RC.get_assembly(wildcards.progenitor),
        bin_bed="results/healr/input_dir/progenitors/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
    output:
        gc_bed="results/healr/input_dir/progenitors/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_gc.bed",
    log:
        "results/logs/bin_information/gc_content/{progenitor}.log",
    params:
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell:
        """
        bedtools nuc -fi {input.assembly} -bed {input.bin_bed} > {output}
        """


def get_all_fastqc_files(wildcards):

    read_files = ["results/fastqc/" + read_path + "_fastqc.zip" for read_path in all_read_paths]

    return read_files


rule compute_mappability:
    input:
        get_all_fastqc_files,
        assembly=lambda wildcards: EAGLE_RC.get_assembly(wildcards.progenitor),
        bin_bed="results/healr/input_dir/progenitors/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
    output:
        genmap_index_dir=directory("results/genmap/{progenitor}/index_dir")
    log:
        "results/logs/bin_information/genmap/{progenitor}.log",
    params:
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell: # Assumes longest read is the read length.
        """
        fastqc_zip_file=$(find "results/fastqc/" -type f -name "*_fastqc.zip" | head -n 1)
        unzip ${{fastqc_zip_file}} -d results/genmap/tmp_{wildcards.progenitor}/
        read_length=$(find results/genmap/tmp_{wildcards.progenitor}/ -type f -name "fastqc_data.txt" | xargs -I{{}} grep "Sequence length" {{}} | awk '{{print $3}}' | awk -F'-' '{{print $2}}')
        rm -r results/genmap/tmp_{wildcards.progenitor}/

        
        mkdir -p $(dirname {output.genmap_index_dir})
        genmap index -F {input.assembly} -I {output.genmap_index_dir}
        genmap map -K ${{read_length}} -E 0 -I {output.genmap_index_dir} -O {output.genmap_index_dir} -bg
        
        genmap_bedgraph=$(find {output.genmap_index_dir} -name *.genmap.bedgraph)
        mappa_bed="results/healr/input_dir/progenitors/{wildcards.progenitor}/{wildcards.progenitor}_${{read_length}}kmer_{params.bin_size}_mappability.bed"
        
        bedtools intersect -a {input.bin_bed} -b $genmap_bedgraph -wo > results/genmap/tmp_{wildcards.progenitor}.genmap.intersect.bed
        bash {workflow.basedir}/scripts/map_in_bins.sh {params.bin_size} results/genmap/tmp_{wildcards.progenitor}.genmap.intersect.bed results/genmap/tmp_{wildcards.progenitor}.genmap.unsorted.bed
        bedtools sort -i results/genmap/tmp_{wildcards.progenitor}.genmap.unsorted.bed > ${{mappa_bed}}
        rm -r results/genmap/tmp_{wildcards.progenitor}*
        """