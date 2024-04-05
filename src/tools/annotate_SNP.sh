#!/bin/bash
#bedtools are required for the task
print_help(){
echo  "
    This script help you locate SNP and transform it to rsid.
    Only intersected SNPs will be printed out.
    Note!!!! 1. All file path should be the same as the original files. 2. bedtools is required, so make sure it's in your environment or path
    Here are the required parameters:
    -g | GWAS summary statistics files. In case you would like to directly combine the output file with a GWAS summary statistics table 
        (Yeah, the main purpose for the script is to annotate the GWAS table). Or any other files that you would like to annotate SNPs to.
    -o | Outpur directory and file header. If you give the value as dir1/dir2/output, then the final output file will be dir1/dir2/output.annot.table.txt
    -c | Which column should the chromosome be (in the GWAS summary table), if parameter 'g' is specified, this must be a given. Or for default it will be 1
    -p | Which column should the position be (in the GWAS summary table), if parameter 'g' is specified, this must be a given. Or for default it will be 2
    -b | to specify what your gwas file starts with 0 or starts with 1. The default parameters are start with 1. #only these two are available!
    -k | Is it necessary to keep all left SNPs although they are not mapped to the dbSNP file? Default is True. (You can choose True or False)
    
    optional parameters:
    -v | VCF file for the dbsnp version you used. For default it will be located in data/ref/Hg19.dbsnp151.vcf. This can be downloaded online
    "
}
#read in parameters
src_dir=$(dirname "$0")
base_dir=$(cd "$src_dir"/../ || exit; pwd)
vcf_file=${base_dir}/data/ref/Hg19.dbsnp151.vcf
chr_col=1
pos_col=2
base_num=1
keep_all="True"
g="None"
while getopts "o:v:g:c:p:b:k:h" opt; do
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
  esac
done
#create snp loc bed file for the input
if [ ! -f ${gwas_file} ]; then
  echo "GWAS file ${gwas_file} invalid!"
  exit
fi
if [ ! -f ${vcf_file} ]; then
  echo "VCF file ${vcf_file} invalid!"
  exit
fi
if [ ${base_num} -ne 1 ] && [ ${base_num} -ne 2 ]; then
  echo "Error: the parameter -b must be 1 or 2"
  exit 
fi
if [ ${keep_all} != "True" ] && [ ${keep_all} != "False" ]; then
  echo "Error: the parameter -k must be either True or False"
  exit 
fi

ori_num=`wc -l ${gwas_file}`

#make a temporary directory
mkdir -p ${output_file}

#make a bed file
awk -v c=${chr_col} -v p=${pos_col} -v b=${base_num} 'BEGIN {OFS="\t"; print "chrom", "from", "to"} {if (NR>1) print $c, $p-b, $p+1-b }' ${gwas_file} > ${output_file}/snp_loc.bed

#get the intersected snp
bedtools intersect -a ${output_file}/snp_loc.bed  -b ${vcf_file} -wa -wb > ${output_file}/temp_loc.out

#get unique output for snp map 
awk '!seen[$1, $2]++' ${output_file}/temp_loc.out > ${output_file}.temp_loc.uniq.out #unique column of rsid
rm -r ${output_file}
fin_num=`wc -l ${output_file}.temp_loc.uniq.out | cut -d " " -f 1`
echo "The original SNP table has ${ori_num} records. After searching for the vcf file, there are ${fin_num} named after rs*."

#merge with the original table
awk 'NR==1{OFS="\t";print $0, "marker.name"}' ${gwas_file} > ${output_file}.annot.table
vcf_base_loc=$((base_num+2))
awk -v c=${chr_col} -v p=${pos_col} -v v=${vcf_base_loc} -v k=${keep_all} '
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
' ${output_file}.temp_loc.uniq.out ${gwas_file} >> ${output_file}.annot.table
rm ${output_file}.temp_loc.uniq.out

#end
echo "Final output file is already in ${output_file}.annot.table"

