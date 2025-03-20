rule define_bins:
    input:
        assembly=EAGLE-RC.get_assembly(PROGENITORS),
    output:
        bin_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
    log:
        "results/logs/bin_information/define_bins_{progenitor}.log",
    params:
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell:
        """
        samtools faidx {input.assembly} 
        bedtools makewindows -g {input.assembly}.fai -w {params.bin_size} > {output.bin_bed}
        """

rule count_gc:
    input:
        assembly=EAGLE-RC.get_assembly(PROGENITORS),
        bin_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
    output:
        gc_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_gc.bed",
    log:
        "results/logs/bin_information/gc_content_{progenitor}.log",
    params:
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell:
        """
        bedtools nuc -fi {input.assembly} -bed {input.bin_bed} > {output}
        """


rule compute_mappability:
    input:
        assembly=EAGLE-RC.get_assembly(PROGENITORS),
        bin_bed="results/bin_information/{progenitor}/{progenitor}_"+f"{BIN_SIZE}_bins.bed",
        fastqc_out_dir="results/fastqc"
    output:
        genmap_index_dir="results/genmap/index_dir"
    log:
        "results/logs/bin_information/genmap_{progenitor}.log",
    params:
        genmap_indexed_assemblies=get_genmap_indexed_assemblies(PROGENITORS),
        bin_size=f"{BIN_SIZE}",
    conda:
        "../envs/bins_gc_map.yaml"
    shell:
        """
        fastqc_zip_dir=$(find {input.fastqc_out_dir} -type f -name "*_fastqc.zip" | head -n 1)
        unzip $fastqc_zip_dir -d results/genmap/tmp
        read_length=$(find results/genmap/tmp -type f -name "fastqc_data.txt" | grep "Sequence length" | awk '{print $3}')
        rm -r results/genmap/tmp

        
        genmap index -F {input.assembly} -I {output.genmap_index_dir}
        genmap map -K $read_length -E 0 -I {output.genmap_index_dir} -O {output.genmap_index_dir} -bg
        
        genmap_bedgraph=$(find {output.genmap_index_dir} -name *.genmap.bedgraph)
        mappa_bed="results/bin_information/{progenitor}/{progenitor}_${read_length}kmer_mappability.bed"
        
        bedtools intersect -a {input.bin_bed} -b $genmap_bedgraph -wo > results/genmap/tmp.genmap.intersect.bed
        bash {workflow.basedir}/scripts/map_in_bins.sh {params.bin_size} results/genmap/tmp.genmap.intersect.bed results/genmap/.tmp.genmap.unsorted.bed
        bedtools sort -i results/genmap/.tmp.genmap.unsorted.bed > ${mappa_bed}
        rm results/genmap/.tmp/genmap.*
        """