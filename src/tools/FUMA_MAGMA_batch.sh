#!/bin/bash
#this file was used to perform fuma and magma analysis
print_help(){
echo  "
    All path should be relative path to the base directory (aka the parent directory of the bash script file).
    Here are the required parameters:
    -m | Mode: FUMA ('fuma') or MAGMA ('magma') analysis
    -e | Expression file.
    -g | GWAS directory. Processed GWAS files are required (*.genes.raw). Where * will be used as the hader
    -o | output file directory
    
    optional parameters:
    -p | MAGMA software directory. The default will be bin/magma/magma.
    -n | header name
    "
}
base_dir=$(dirname "$(dirname "$(dirname "$(realpath "$0")")")")
magma_path=${base_dir}/bin/magma/magma

while getopts "m:e:g:p:o:n:h" opt; 
do
  case $opt in
    h)
      print_help
      exit 0 ;;
    m)
      mode=$OPTARG ;;
    e)
      exp_file=${base_dir}/$OPTARG ;;
    g)
      gwas_dir=${base_dir}/$OPTARG ;;
    o)
      output_dir=${base_dir}/$OPTARG ;;
    p)
      magma_path=${base_dir}/$OPTARG ;;
    #n)
    #  header_name=$OPTARG ;;
    *)
      print_help
      exit 0 ;;
  esac
done

measure_time() {
  local cmd="$1"
  local start="$(date +%s%3N)" # start time in milliseconds
  eval "${cmd}"
  local end=$(date +%s%3N)   # end time in milliseconds
  local duration_ms=$((end-start))
  local duration_s=$(echo "scale=3; ${duration_ms} / 1000" | bc)
  echo "Running time for \"${cmd}\": ${duration_s} s"
}

#check
if [ ! -f "${exp_file}" ]; then
  echo "expression file ${exp_file} invalid!"
  exit
fi
if [ ! -d "${gwas_dir}" ] || ! echo "${gwas_dir}"/* | grep -q ".genes.raw"; then
  echo "The gwas_directory doesn't exists or there aren't any genes.raw file there！"
  exit
fi

#make up output directory
if [ ! -d "${output_dir}" ]; then
  mkdir -p "${output_dir}"
fi

#time
time_string=()
#magma or fuma
if [ "$mode" = "MAGMA" ] || [ "$mode" = "Magma" ] || [ "$mode" = "magma" ]; then
  for gwas_file in "$gwas_dir"/*.genes.raw; do
    base_name=$(basename "$gwas_file")
    base_name=${base_name%.genes.raw}
    start=$(date +%s%3N)
    nohup $magma_path --gene-results "$gwas_file" --set-annot "$exp_file" --out "$output_dir"/"$base_name" &
    end=$(date +%s%3N)
    duration=$((end-start))
    time_string+=("for ${gwas_file} is ${duration} ms")
  done
elif [ "$mode" = "FUMA" ] || [ "$mode" = "Fuma" ] || [ "$mode" = "fuma" ]; then
  for gwas_file in "$gwas_dir"/*.genes.raw; do
    base_name=$(basename "$gwas_file")
    base_name=${base_name%.genes.raw}
    start=$(date +%s%3N)
    nohup $magma_path --gene-results "$gwas_file" --gene-covar "$exp_file" --model condition-hide=Average direction=greater --out "$output_dir"/"$base_name" &
    end=$(date +%s%3N)
    duration=$((end-start))
    time_string+=("for ${gwas_file} is ${duration} ms")
  done
else
  echo "Error, the mode should be MAGMA or FUMA, you are specifying ${mode}"
  exit
fi

for i in "${time_string[@]}";do
  echo "$i"
done