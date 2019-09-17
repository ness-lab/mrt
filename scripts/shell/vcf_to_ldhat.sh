
vcftools --vcf $1 --chr 1 --ldhat &&

filename=$(basename -- "${1%.*}")

rm out.log

mkdir -p ldhat_input_files/

mv out.ldhat.sites "ldhat_input_files/${filename}.sites"
mv out.ldhat.locs "ldhat_input_files/${filename}.locs"

