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
        samtools faidx {input} 2> {log}
        bedtools makewindows -g {input}.fai -w {params.bin_size} > {output.bin_bed} 2> {log}
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
        bedtools nuc -fi {input.assembly} -bed {input.bin_bed} > {output} 2> {log}
        """

def get_random_read_file(filter):
    
    if filter == True:
        polyploids=glob.glob(f"{INPUT_DIR}/polyploids/*")
        random_sample_name=os.path.basename(polyploids[0])
        random_file = f"results/fastp/{random_sample_name}/{random_sample_name}_R1_filtered.fastq"
        return random_file
    else:    
        return os.path.join(f"{INPUT_DIR}/polyploids/", all_read_paths[0]),


rule compute_mappability:
    input:
        random_read_file=get_random_read_file(FILTER),
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
    shell:# Assumes all short read files have same length of reads
        """

        random_read_file="{input.random_read_file}"
        read_length=$(
            {{
                seqkit head -n 500 "${{random_read_file}}" \
                | seqkit fx2tab -l \
                | awk '{{print $4}}' \
                | sort | uniq -c | sort -nr | head -1 \
                | awk '{{print $2}}' 
            }} 2>&1 | tee -a "{log}"
        )
        
        mkdir -p "results/genmap/tmp_{wildcards.progenitor}" 2>&1 | tee -a "{log}" 
        mkdir -p "$(dirname "{output.genmap_index_dir}")" 2>&1 | tee -a "{log}"
        genmap index -F "{input.assembly}" -I "{output.genmap_index_dir}" 2>&1 | tee -a "{log}"
        genmap map -K "${{read_length}}" -E 0 -I "{output.genmap_index_dir}" -O "{output.genmap_index_dir}" -bg 2>&1 | tee -a "{log}"
        
        genmap_bedgraph=$(find "{output.genmap_index_dir}" -name *.genmap.bedgraph)
        mappa_bed="results/healr/input_dir/progenitors/{wildcards.progenitor}/{wildcards.progenitor}_${{read_length}}kmer_{params.bin_size}_mappability.bed"
        
        bedtools intersect -a "{input.bin_bed}" -b "${{genmap_bedgraph}}" -wo > "results/genmap/tmp_{wildcards.progenitor}.genmap.intersect.bed" 2>&1 | tee -a "{log}"
        bash "{workflow.basedir}/scripts/map_in_bins.sh" "{params.bin_size}" "results/genmap/tmp_{wildcards.progenitor}.genmap.intersect.bed" "results/genmap/tmp_{wildcards.progenitor}.genmap.unsorted.bed" mkdir -p "$(dirname "{output.genmap_index_dir}")" 2>&1 | tee -a "{log}"
        bedtools sort -i "results/genmap/tmp_{wildcards.progenitor}.genmap.unsorted.bed" > "${{mappa_bed}}" 2>&1 | tee -a "{log}"
        rm -r "results/genmap/tmp_{wildcards.progenitor}"* 2>&1 | tee -a "{log}"
        """