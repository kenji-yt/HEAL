#########################################
## Let's make a reproducibility report ##
#########################################

input_dir=$1
n_cores=$2
snake_genespace_report=results/snake_GENESPACE_reproducibility_report.txt
snake_eagle_rc_report=results/snake_EAGLE_RC_reproducibility_report.txt
report=results/HEAL_reproducibility_report.txt
CURRENT_DATETIME=$(date +"%Y-%m-%d %H:%M:%S")

echo "**************" >> "${report}"
echo "*    HEAL    *" >> "${report}"
echo "**************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "Reproducibility report for HEAL." >> "${report}"
echo "Run date & time: ${CURRENT_DATETIME}" >> "${report}"
echo "Number of allocated cores: ${n_cores}" >> "${report}" 
echo "" >> "${report}"
echo "" >> "${report}"


version_HEAL=$(git describe --tags --abbrev=0 | sed 's/v//g')
echo "HEAL=${version_HEAL}" >> "${report}"
