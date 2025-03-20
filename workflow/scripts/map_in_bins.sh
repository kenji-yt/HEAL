bin_size=$1
intersect_bed=$2
unsorted_bed=$3

awk -v bin_size=${bin_size} '{
            bin_start = $2
            bin_end = $3
            overlap_length = $NF
            mappability = $(NF-1)
            
            bin_key = $1 "_" bin_start "_" bin_end
            weighted_sum[bin_key] += mappability * overlap_length
        }
        END {
            for (bin in weighted_sum) {
                split(bin, bin_coords, "_")
                print bin_coords[1]"\t"bin_coords[2]"\t"bin_coords[3]"\t"weighted_sum[bin] / bin_size
            }
        }' ${intersect_bed} > ${unsorted_bed}
