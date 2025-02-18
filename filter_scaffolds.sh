#!/bin/bash

#get input:

#no paths used, make sure to run in the directory where all the assembly output follders are stored. phi6 refs have to be stored here too for blast

#-p ${parent_dir} directory containing barcode folders
#-d sample_suffix - starts with an underscore, example: "_norm3_meta_spades"
#-s ${host_s} S segment reference file & path
#-m ${host_m} M segment reference file & path
#-l ${host_l} L segment reference file & path


#parent_path=''
sample_suffix=''
#host_s=''
#host_m=''
#host_l=''


print_usage() {
  printf "Usage: ..."
}

while getopts d: flag
do
    case "${flag}" in
d) sample_suffix=${OPTARG};;
#s) host_s=${OPTARG};;
#m) host_m=${OPTARG};;
#l) host_l=${OPTARG};;

    esac
done


#cd ${parent_path} || exit

echo "working shell directory:"
pwd
echo


# #Blast the S segment 
# for direct in *${sample_suffix}_output; 
# do file=$(basename $direct); 
# sample=${file##*/}
# sample=${file%%_*}; 
# #can add line here to remove contig header whitespace if necessary
# #blast the current assembly against phi6 M segment reference
# blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_S.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv;
#         blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv)"
#         if [ -n "${blast_check}" ]; then
#                 echo "${sample} ${blast_check} not empty"
#                 #makes file with scaffold names only
#                 awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout_names.txt;
#                 cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_scaffolds.fasta
#         fi;
# done


# #Blast the M segment 
# for direct in *${sample_suffix}_output; 
# do file=$(basename $direct); 
# sample=${file##*/}
# sample=${file%%_*}; 
# #can add line here to remove contig header whitespace if necessary
# #blast the current assembly against phi6 M segment reference
# blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_M.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv;
#         blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv)"
#         if [ -n "${blast_check}" ]; then
#                 echo "${sample} ${blast_check} not empty"
#                 #makes file with scaffold names only
#                 awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout_names.txt;
#                 cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_scaffolds.fasta
#         fi;
# done

# #Blast the L segment 
# for direct in *${sample_suffix}_output; 
# do file=$(basename $direct); 
# sample=${file##*/}
# sample=${file%%_*}; 
# #can add line here to remove contig header whitespace if necessary
# #blast the current assembly against phi6 M segment reference
# blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_L.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv;
#         blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv)"
#         if [ -n "${blast_check}" ]; then
#                 echo "${sample} ${blast_check} not empty"
#                 #makes file with scaffold names only
#                 awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout_names.txt;
#                 cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_scaffolds.fasta
#         fi;
# done


#S segment matches filtering!
#print output to file with the following headers: 
echo "sample,segment,num scaffolds,max scaffold length,min scaffold length,max kmer coverage,min kmer coverage" > S_all_scaffold_stats.csv

for file in *${sample_suffix}_output/*${sample_suffix}_S_scaffolds.fasta; 
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
    cov_list+=("${cov}");    

    #extract length value (as string)
    len=${scaff##*length_};
    len=${len%%_*};
    #add kmer ov value to list
    len_list+=("${len}");
    

  done;
  #sort the list based on numerical value, largest to smallest
  sorted_cov_list=(`printf '%s\n' "${cov_list[@]}" | sort -gr`); 
  #sort the list based on numerical value, largest to smallest
  sorted_len_list=(`printf '%s\n' "${len_list[@]}" | sort -gr`); 
  #echo "len values = ${sorted_len_list[@]} for $sample";
  
  #calculate a % of the max cov value:
  thirty_perc_cov=$(echo "${sorted_cov_list[0]} * 0.3" | bc);
  echo "30% of max cov = $thirty_perc_cov";
  filtered_cov_values=();
  sort_filt_cov_values=();

  #calculate %of the max length value:
  twenty_perc_len=$(echo "${sorted_len_list[0]} * 0.2" | bc);
  echo "20% of max len = $twenty_perc_len";
  filtered_len_values=();
  sort_filt_len_values=();

  for i in $(seq 0 $( echo "$(grep -c '>' $file)" | bc)); 
    do if (( $(echo "${cov_list[i]} > 15" |bc -l) )); then
      echo "${len_list[i]} & ${cov_list[i]}";
      if ( (( $(echo "${cov_list[i]} > $thirty_perc_cov" |bc -l) )) || (( $(echo "${cov_list[i]} > 300" |bc -l) )) ) && ( (( $(echo "${len_list[i]} > $twenty_perc_len" |bc -l) )) || (( $(echo "${len_list[i]} > 800" |bc -l) )) ); then
        filtered_cov_values+=("${cov_list[i]}");
        filtered_len_values+=("${len_list[i]}");
        grep ${cov_list[i]} $file >> ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt;
      fi;
    fi;
  done;
  sort_filt_cov_values=(`printf '%s\n' "${filtered_cov_values[@]}" | sort -gr`); 
  sort_filt_len_values=(`printf '%s\n' "${filtered_len_values[@]}" | sort -gr`);
  #after going thru whole file, add the sequences baseed on the filtered headers to new fasta called ..covfilt...fasta
  awk 'BEGIN{FS="\n";RS=">";ORS=""} (NR==FNR)&&(NR>1){headers[$1];next} ($1 in headers){print ">" $0}' ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds.fasta;
  echo "kcov values above 30% or 300 = ${filtered_cov_values[@]} + length values above 20% or 800 = ${filtered_len_values[@]}";
  echo;
  echo "${sample},S,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> S_all_scaffold_stats.csv;
done



#filter M segment matches!
#print output to file with the following headers: 
echo "sample,segment,num scaffolds,max scaffold length,min scaffold length,max kmer coverage,min kmer coverage" > M_all_scaffold_stats.csv

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
    cov_list+=("${cov}");
     
    #extract length value (as string)
    len=${scaff##*length_};
    len=${len%%_*};
    #add kmer ov value to list
    len_list+=("${len}");
  done;

  #sort the list based on numerical value, largest to smallest
  sorted_cov_list=(`printf '%s\n' "${cov_list[@]}" | sort -gr`); 
  #sort the list based on numerical value, largest to smallest
  sorted_len_list=(`printf '%s\n' "${len_list[@]}" | sort -gr`); 
  #echo "len values = ${sorted_len_list[@]} for $sample";
    
  #calculate a % of the max cov value:
  thirty_perc_cov=$(echo "${sorted_cov_list[0]} * 0.3" | bc);
  echo "30% of max cov = $thirty_perc_cov";
  filtered_cov_values=();
  sort_filt_cov_values=();

  #calculate %of the max length value:
  twenty_perc_len=$(echo "${sorted_len_list[0]} * 0.2" | bc);
  echo "20% of max len = $twenty_perc_len";
  filtered_len_values=();
  sort_filt_len_values=();

  for i in $(seq 0 $( echo "$(grep -c '>' $file)" | bc)); 
    do if (( $(echo "${cov_list[i]} > 15" |bc -l) )); then
      echo "${len_list[i]} & ${cov_list[i]}";
      if ( (( $(echo "${cov_list[i]} > $thirty_perc_cov" |bc -l) )) || (( $(echo "${cov_list[i]} > 300" |bc -l) )) ) && ( (( $(echo "${len_list[i]} > $twenty_perc_len" |bc -l) )) || (( $(echo "${len_list[i]} > 800" |bc -l) )) ); then
        filtered_cov_values+=("${cov_list[i]}");
        filtered_len_values+=("${len_list[i]}");
        grep ${cov_list[i]} $file >> ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt;
      fi;
    fi;
  done;
  sort_filt_cov_values=(`printf '%s\n' "${filtered_cov_values[@]}" | sort -gr`); 
  sort_filt_len_values=(`printf '%s\n' "${filtered_len_values[@]}" | sort -gr`);
  #after going thru whole file, add the sequences baseed on the filtered headers to new fasta called ..covfilt...fasta
  awk 'BEGIN{FS="\n";RS=">";ORS=""} (NR==FNR)&&(NR>1){headers[$1];next} ($1 in headers){print ">" $0}' ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds.fasta;
  echo "kcov values above 30% or 300 = ${filtered_cov_values[@]} + length values above 20% or 800 = ${filtered_len_values[@]}";
  echo;
  echo "${sample},M,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> M_all_scaffold_stats.csv;
done



#L segment matches filtering!
#print output to file with the following headers: 
echo "sample,segment,num scaffolds,max scaffold length,min scaffold length,max kmer coverage,min kmer coverage" > L_all_scaffold_stats.csv

for file in *${sample_suffix}_output/*${sample_suffix}_L_scaffolds.fasta; 
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
    cov_list+=("${cov}");
    #extract length value (as string)
    len=${scaff##*length_};
    len=${len%%_*};
    #add kmer ov value to list
    len_list+=("${len}");
    

  done;
  #sort the list based on numerical value, largest to smallest
  sorted_cov_list=(`printf '%s\n' "${cov_list[@]}" | sort -gr`); 
  #sort the list based on numerical value, largest to smallest
  sorted_len_list=(`printf '%s\n' "${len_list[@]}" | sort -gr`); 
  #echo "len values = ${sorted_len_list[@]} for $sample";
  
  #calculate a % of the max cov value:
  thirty_perc_cov=$(echo "${sorted_cov_list[0]} * 0.3" | bc);
  echo "30% of max cov = $thirty_perc_cov";
  filtered_cov_values=();
  sort_filt_cov_values=();

  #calculate %of the max length value:
  twenty_perc_len=$(echo "${sorted_len_list[0]} * 0.2" | bc);
  echo "20% of max len = $twenty_perc_len";
  filtered_len_values=();
  sort_filt_len_values=();

  for i in $(seq 0 $( echo "$(grep -c '>' $file)" | bc)); 
    do if (( $(echo "${cov_list[i]} > 15" |bc -l) )); then
      echo "${len_list[i]} & ${cov_list[i]}"; 
      if ( (( $(echo "${cov_list[i]} > $thirty_perc_cov" |bc -l) )) || (( $(echo "${cov_list[i]} > 300" |bc -l) )) ) && ( (( $(echo "${len_list[i]} > $twenty_perc_len" |bc -l) )) || (( $(echo "${len_list[i]} > 800" |bc -l) )) ); then
        filtered_cov_values+=("${cov_list[i]}");
        filtered_len_values+=("${len_list[i]}");
        grep ${cov_list[i]} $file >> ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt;
      fi;
    fi;
  done;
  sort_filt_cov_values=(`printf '%s\n' "${filtered_cov_values[@]}" | sort -gr`); 
  sort_filt_len_values=(`printf '%s\n' "${filtered_len_values[@]}" | sort -gr`);
  #after going thru whole file, add the sequences baseed on the filtered headers to new fasta called ..covfilt...fasta
  awk 'BEGIN{FS="\n";RS=">";ORS=""} (NR==FNR)&&(NR>1){headers[$1];next} ($1 in headers){print ">" $0}' ${sample}${sample_suffix}_output/${fileonly}_covlenfilt_names.txt ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds.fasta;
  echo "kcov values above 30% or 300 = ${filtered_cov_values[@]} + length values above 20% or 800 = ${filtered_len_values[@]}";
  echo;
  echo "${sample},L,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> L_all_scaffold_stats.csv;
done
