

# Make new directory to store sorted VCFs
new_dir=${1}/sorted_vcfs
mkdir -p $new_dir

# Loop through VCFs, sort, and add to new directory
for file in $1/*.vcf
do
  filename=$(basename -- "${file%.*}")
  extension="${file##*.}"
  new_name="${filename}_sorted.${extension}"

  # Sort by second filed in VCF (i.e. POS)
  sort -n -k 2 $file > ${new_dir}/${new_name}
done
