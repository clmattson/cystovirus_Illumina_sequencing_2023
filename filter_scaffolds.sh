#!/bin/bash

#get input:

#-p ${parent_dir} directory containing barcode folders
#-d sample_suffix - starts with an underscore
#-s ${host_s} S segment reference file & path
#-m ${host_m} M segment reference file & path
#-l ${host_l} L segment reference file & path


parent_path=''
sample_suffix=''
host_s=''
host_m=''
host_l=''


print_usage() {
  printf "Usage: ..."
}

while getopts p:d:s:m:l: flag
do
    case "${flag}" in
p) parent_path=${OPTARG};;
d) sample_suffix=${OPTARG};;
s) host_s=${OPTARG};;
m) host_m=${OPTARG};;
l) host_l=${OPTARG};;

    esac
done


cd ${parent_path} || exit

echo "working shell directory:"
pwd
echo




for direct in *${sample_suffix}_output; 
do file=$(basename $direct); 
sample=${file##*/}
sample=${file%%_*}; 
#can add line here to remove contig header whitespace if necessary
#blast the current assembly against phi6 M segment reference
blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_M.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds_M_blastout.csv;
        blast_check="$(cat ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds_M_blastout.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${sample} ${blast_check} not empty"
                #makes file with scaffold names only
                awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds_M_blastout.csv > ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds_M_blastout_names.txt;
                cat ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds_M_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}_norm3_meta_spades_M_scaffolds.fasta
        fi;
done





#print output to file with the following headers: 
echo "sample,num scaffolds,max scaffold length,min scaffold length,max kmer coverage,min kmer coverage" > M_all_scaffold_stats.csv

for file in *${sample_suffix}_output/*${sample_suffix}_M_scaffolds.fasta; 
  do cov_list=();
  len_list=();
  #get filename
  fileonly=$(basename $file); 
  #get sample letter from file name:
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
  echo "finding scaffolds with kmer cov > 30% of max for sample ${sample}";
  #loop thru the scaffold headers in the current sample scaffold file:
  for scaff in $(grep '>' $file); 
    #extract coverage value (as string)
    do cov=${scaff##*_};
    echo "current scaffold = ${scaff}"
    #add kmer ov value to list
    if (( $(echo "$cov > 15" |bc -l) )); then
      cov_list+=("${cov}");
    fi;
    #sort the list based on numerical value, largest to smallest
    sorted_cov_list=(`printf '%s\n' "${cov_list[@]}" | sort -gr`); 

    #extract length value (as string)
    len=${scaff##*length_};
    len=${len%%_*};
    #add kmer ov value to list
    len_list+=("${len}");
    #sort the list based on numerical value, largest to smallest
    sorted_len_list=(`printf '%s\n' "${len_list[@]}" | sort -gr`); 
    #echo "len values = ${sorted_len_list[@]} for $sample";

  done;
  #calculate a % of the max cov value:
  thirty_perc_cov=$(echo "${sorted_cov_list[0]} * 0.3" | bc);
  echo "30% of max cov = $thirty_perc_cov";
  filtered_cov_values=();

  #calculate %of the max length value:
  twenty_perc_len=$(echo "${sorted_len_list[0]} * 0.2" | bc);
  echo "20% of max len = $twenty_perc_len";
  filtered_len_values=();

  for i in $(seq 0 $( echo "$(grep -c '>' $file)" | bc)); 
    do if ( (( $(echo "${cov_list[i]} > $thirty_perc_cov" |bc -l) )) || (( $(echo "${cov_list[i]} > 300" |bc -l) )) ) && ( (( $(echo "${len_list[i]} > $twenty_perc_len" |bc -l) )) || (( $(echo "${len_list[i]} > 800" |bc -l) )) ); then
      filtered_cov_values+=("${cov_list[i]}");
      filtered_len_values+=("${len_list[i]}");
      grep ${cov_list[i]} $file >> *${sample_suffix}_output/${fileonly}_covlenfilt_names.txt;
    fi;
  done;
  #after going thru whole file, add the sequences baseed on the filtered headers to new fasta called ..covfilt...fasta
  awk 'BEGIN{FS="\n";RS=">";ORS=""} (NR==FNR)&&(NR>1){headers[$1];next} ($1 in headers){print ">" $0}' *${sample_suffix}_output/${fileonly}_covlenfilt_names.txt *${sample_suffix}_output/${sample}_norm3_meta_spades_M_scaffolds.fasta >> *${sample_suffix}_output/${sample}_norm3_meta_spades_M_covlenfilt_scaffolds.fasta;
  echo "kcov values above 30% or 300 = ${filtered_cov_values[@]} + length values above 20% or 800 = ${filtered_len_values[@]}";
  echo;
  echo "${sample},$(grep -c '>' *${sample_suffix}_output/${sample}_norm3_meta_spades_M_covlenfilt_scaffolds.fasta),${filtered_len_values[0]},${filtered_len_values[-1]},${filtered_cov_values[0]},${filtered_cov_values[-1]}" >> M_all_scaffold_stats.csv;
done
