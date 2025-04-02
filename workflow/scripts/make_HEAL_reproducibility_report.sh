#########################################
## Let's make a reproducibility report ##
#########################################

input_dir=$1
n_cores=$2
snake_eagle_rc_report=results/snake_EAGLE_RC_reproducibility_report.txt
snake_genespace_report=results/snake_GENESPACE_reproducibility_report.txt
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


echo "********************" >> "${report}"
echo "* Operating System *" >> "${report}"
echo "********************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

OS=$(uname -s)

if [ "$OS" == "Linux" ]; then
    # For Linux, try to get version from /etc/os-release
    if [ -f /etc/os-release ]; then
        source /etc/os-release
        echo "Operating System: $NAME" >> "${report}"
        echo "Version: $VERSION" >> "${report}"
    else
        echo "Linux OS (version unknown)" >> "${report}"
    fi
# Assume anything else is macOS
else
    echo "Operating System: macOS"  >> "${report}"
    sw_vers  >> "${report}"
fi

echo "" >> "${report}"
echo "" >> "${report}"
echo "**************" >> "${report}"
echo "* INPUT DATA *" >> "${report}"
echo "**************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"


awk '/\* INPUT DATA \*/ {flag=1; next} /\* TOOLS \*/ {flag=0} flag' $snake_eagle_rc_report \
 | grep -v '^$' | grep -v '^*' >> "${report}"

awk '/\* INPUT DATA \*/ {flag=1; next} /\* TOOLS \*/ {flag=0} flag' $snake_genespace_report \
 | grep -v '^$' | grep -v '^*' >> "${report}"

echo "" >> "${report}"
echo "" >> "${report}"
echo "*********" >> "${report}"
echo "* TOOLS *" >> "${report}"
echo "*********" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

version_HEAL=$(git describe --tags --abbrev=0 | sed 's/v//g')
echo "HEAL=${version_HEAL}" >> "${report}"

awk '/\* TOOLS \*/ {flag=1; next} /\* OUTPUT FILES \*/ {flag=0} flag' $snake_eagle_rc_report \
 | grep -v '^$' | grep -v '^*' >> "${report}"

awk '/\* TOOLS \*/ {flag=1; next} /\* OUTPUT FILES \*/ {flag=0} flag' $snake_genespace_report \
 | grep -v '^$' | grep -v '^*' >> "${report}"

echo "" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"
echo "****************" >> "${report}"
echo "* OUTPUT FILES *" >> "${report}"
echo "****************" >> "${report}"
echo "" >> "${report}"
echo "" >> "${report}"

if [ "$OS" == "Linux" ]; then
    echo "Linux md5sum checksums for the output files" >> "${report}"
    find results/healr/stats/ -type f | xargs -n${n_cores} md5sum | awk '{print $2"\t"$1}' >> "${report}"  
else
    echo "Mac md5 checksums for the input files" >> "${report}"
    find results/healr/stats/ -type f | xargs -n${n_cores} md5 | awk '{print $2"\t"$4}' >> "${report}"  
fi
echo "" >> "${report}"
awk '/\* OUTPUT FILES \*/ {flag=1; next} flag' $snake_eagle_rc_report \
| grep -v '^$' | grep -v '^*' | grep -v '^Mac md5' | grep -v '^Linux md5sum' >> "${report}"
echo "" >> "${report}"
awk '/\* OUTPUT FILES \*/ {flag=1; next} flag' $snake_genespace_report \
| grep -v '^$' | grep -v '^*' | grep -v '^Mac md5' | grep -v '^Linux md5sum' >> "${report}"

rm results/snake_*_reproducibility_report.txt