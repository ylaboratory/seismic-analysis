#!/bin/bash
#bedtools are required for the task
print_help(){
echo  "
    This script help you locate SNP and transform the locus (by chromosome id and location) to the corresponding variant rsid.
    Only intersected SNPs will be printed out.
    Note!!!! 1. All file path should be the same as the original files. 2. bedtools is required, so make sure it's in your environment or path
    Here are the required parameters:
    -g | GWAS summary statistics files. Or any other files that you would like to annotate SNPs to.
    -o | Outpur file header. If the value is given as 'dir1/dir2/output', then the final output file will be 'dir1/dir2/output.annot.table.txt'. If it is not specified, the output file will be the same as the input file with the suffix ignored.
    -c | Chromosome column, i.e. the column with with chromosome information be (in the table). The value should be specified or by default it will be 1.
    -p | Position column, i.e. The column with the the position information (in the table), The value should be specified or by default it will be 2.
    -b | Basis of chromosome numbering, to specify the if first position on a chromosome is labeled as '1'(1-base), or '0'(0-base). The default parameters are start with 1. #only these two are available!
    -k | Keep all unmapped SNPs which do not exist in the dbSNP file or not. By default it is True. (You can choose between 'True' or 'False')
    
    optional parameters:
    -v | VCF file for the dbsnp version you used. For default it will be located in 'data/ref/Hg19.dbsnp151.vcf'. This can be downloaded online.
    -n | Number of tasks to split the job. The default value is 1. If the memory is not enough, you can split the job into several parts but it will take more time.
    -t | Customized bedtools command. By default it is 'bedtools'. If you do not have it in your environment, you can specify the path to the bedtools command.
    -w | The working directory. By default it is the current directory. If you want to specify the working directory, you can use this parameter.
    -h | Help message.
    "
}
#read in parameters
base_dir=$(dirname "$0")
vcf_file=${base_dir}/data/ref/Hg19.dbsnp151.vcf
chr_col=1
pos_col=2
base_num=1
keep_all="True"
gwas_file="None"
output_file="None"
bedtools_command="bedtools"


while getopts "o:v:g:c:p:b:k:n:h:t" opt; do
  case $opt in
  h)
    print_help
    exit 0
    ;;
  o)
    output_file=${base_dir}/$OPTARG
    ;;
  v)
    vcf_file=${base_dir}/$OPTARG
    ;;
  g)
    gwas_file=${base_dir}/$OPTARG
    ;;
  c)
    chr_col=$OPTARG
    ;;
  p)
    pos_col=$OPTARG
    ;;
  b)
    base_num=$OPTARG
    ;;
  k)
    keep_all=$OPTARG
    ;;
  n)
    num_task=$OPTARG
    ;;
  t)
    bedtools_command=$OPTARG
    ;;
  *)
    echo "Invalid parameters!"
    print_help
    exit 1
    ;;
  esac
done

#parameter check
if [ ! -f "$gwas_file" ]; then
  echo "GWAS file ${gwas_file} invalid! Use the -h parameter to see the help message."
  exit
fi
#if output file is not specified, then use the same name as the input file
if [ "$output_file" == "None" ]; then
  output_file=$(dirname "$gwas_file")/$(basename "$gwas_file" ｜rev|cut -d '.' -f 2-|rev)
fi

#other parameter format check
if [ ! -f "$vcf_file" ]; then
  echo "VCF file ${vcf_file} invalid! Use the -h parameter to see the help message."
  exit
fi

if [ "$base_num" -ne 1 ] && [ "$base_num" -ne 2 ]; then
  echo "Error: the parameter -b must be '1' or '2'. Use the -h parameter to see the help message."
  exit 
fi

if [ "$keep_all" != "True" ] && [ "$keep_all" != "False" ]; then
  echo "Error: the parameter -k must be either 'True' or 'False'. Use the -h parameter to see the help message."
  exit 
fi

ori_num=$(wc -l "$gwas_file")
#make a temporary directory
mkdir -p "$output_file"

#if n > 1 then split the vcf file
if [ "$num_task" -gt 1 ]; then
  vcf_file_dir=$(dirname "$vcf_file")
  split -n l/"$num_task" --numeric-suffixes "$vcf_file"  "$vcf_file_dir"/vcf_split.

  for tmp_vcf_file in  "$vcf_file_dir"/vcf_split.*; do
    vcf_suffix=$(echo "$tmp_vcf_file" | rev | cut -d '.' -f 1 | rev)
    #get the first and last line of the vcf file

    first_line=$(awk '!/^#/{print $0; exit}' "$tmp_vcf_file")
    final_line=$(awk '!/^#/{line=$0} END{print line}' "$tmp_vcf_file")

    first_chr=$(echo "$first_line" | cut -f 1)
    first_pos=$(echo "$first_line" | cut -f 2)
    last_chr=$(echo "$final_line" | cut -f 1)
    last_pos=$(echo "$final_line" | cut -f 2)

    #get the intersected snp for the gwas file based on the vcf file
    awk -v c="$chr_col" -v p="$pos_col" -v b="$base_num" -v fc="$first_chr" -v lc="$last_chr" -v fp="$first_pos" -v lp="$last_pos"  '
      function chr_order(chr) {
        if (chr ~ /^[0-9]+$/) return chr + 0
        else if (chr == "X") return 23
        else if (chr == "Y") return 24
        else if (chr == "MT") return 25
        else return 26
      }
      BEGIN {
      first_chr = chr_order(fc)
      last_chr = chr_order(lc)
      }

      NR==1 || /^#/ {print $0; next}
      {
        curr_chr = chr_order($c)
        if ((curr_chr > first_chr || (curr_chr == first_chr && $p >= fp)) && 
            (curr_chr < last_chr || (curr_chr == last_chr && $p <= lp)))
            print $0
        }
    ' "$gwas_file" > "$output_file"/temp_gwas_"$vcf_suffix".txt
  done

fi

process_part() {
    local gwas_part=$1
    local vcf_part=$2
    local output_part=$3
    local chr_col=$4
    local pos_col=$5
    local base_num=$6
    local keep_all=$7

    # Make a bed file
    awk -v c="$chr_col" -v p="$pos_col" -v b="$base_num" 'BEGIN {OFS="\t"; print "chrom", "from", "to"} {if (NR>1) print $c, $p-b, $p+1-b }' "$gwas_part" > "${output_part}.snp_loc.bed"

    # Get the intersected snp
    "$bedtools_command" intersect -a "${output_part}.snp_loc.bed" -b "$vcf_part" -wa -wb > "${output_part}.temp_loc.out"

    # Get unique output for snp map
    awk '!seen[$1, $2]++' "${output_part}.temp_loc.out" > "${output_part}.temp_loc.uniq.out"

    # Merge with the original table
    awk 'NR==1{OFS="\t";print $0, "marker.name"}' "$gwas_part" > "${output_part}.annot.table"
    vcf_base_loc=$((base_num+2))
    awk -v c="$chr_col" -v p="$pos_col" -v v="$vcf_base_loc" -v k="$keep_all" '
      NR == FNR{
        f2[$1, $v] = $6
        next
      }
        NR != FNR && FNR>1  {
          if (($c, $p) in f2){
            print $0, f2[$c, $p]
          }else if (k=="True"){
            print $0, $c":"$p
          }
         }
    ' "${output_part}.temp_loc.uniq.out" "$gwas_part" >> "${output_part}.annot.table"

    # Clean up temporary files
    rm "${output_part}.snp_loc.bed" "${output_part}.temp_loc.out" "${output_part}.temp_loc.uniq.out"
}

#final pipeline
if [ "$num_task" -gt 1 ]; then
  for tmp_vcf_file in  "$vcf_file_dir"/vcf_split.*; do
    vcf_suffix=$(echo "$tmp_vcf_file" | rev | cut -d '.' -f 1 | rev)
    mkdir -p "$output_file"/temp_gwas_"$vcf_suffix"
    process_part "$output_file"/temp_gwas_"$vcf_suffix".txt "$tmp_vcf_file" "$output_file"/temp_gwas_"$vcf_suffix" "$chr_col" "$pos_col" "$base_num" "$keep_all"
  done
  #merge all the files
  awk 'NR==1{print $0} FNR>1' "$output_file"/temp_gwas_*.annot.table > "$output_file".annot.table
  rm -r "$output_file"/temp_gwas_*
else
  process_part "$gwas_file" "$vcf_file" "$output_file" "$chr_col" "$pos_col" "$base_num" "$keep_all"
fi

#end
echo "Final output file is already in '$output_file'.annot.table"

