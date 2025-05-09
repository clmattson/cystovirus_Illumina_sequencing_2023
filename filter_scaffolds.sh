#!/bin/bash

#to do:
#run script on non-normalized data. maybe switch output to like less long format version?  
#get pariwaise alignment %s for filtered scaffold files and add to info table. 


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


#Blast the S segment 
for direct in *${sample_suffix}_output; 
do file=$(basename $direct); 
sample=${file##*/}
sample=${file%%_*}; 
#can add line here to remove scaffold header whitespace if necessary
#blast the current assembly against phi6 M segment reference
blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_S.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv;
        blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${sample} ${blast_check} not empty"
                #makes file with scaffold names only
                awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout_names.txt;
                cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_S_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_scaffolds.fasta
        fi;
done


#Blast the M segment 
for direct in *${sample_suffix}_output; 
do file=$(basename $direct); 
sample=${file##*/}
sample=${file%%_*}; 
#can add line here to remove scaffold header whitespace if necessary
#blast the current assembly against phi6 M segment reference
blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_M.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv;
        blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${sample} ${blast_check} not empty"
                #makes file with scaffold names only
                awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout_names.txt;
                cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_M_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_scaffolds.fasta
        fi;
done

#Blast the L segment 
for direct in *${sample_suffix}_output; 
do file=$(basename $direct); 
sample=${file##*/}
sample=${file%%_*}; 
#can add line here to remove scaffold header whitespace if necessary
#blast the current assembly against phi6 M segment reference
blastn -query ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta -subject phi6_L.fasta -outfmt '6 delim=,' -qcov_hsp_perc 20 >  ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv;
        blast_check="$(cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${sample} ${blast_check} not empty"
                #makes file with scaffold names only
                awk -F"," '{print $1}' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout.csv > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout_names.txt;
                cat ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds_L_blastout_names.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${sample}${sample_suffix}_output/${sample}${sample_suffix}_scaffolds.fasta >> ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_scaffolds.fasta
        fi;
done


#S segment matches filtering!
#print output to file with the following headers: 
echo "sample,segment,num S scaffolds,max S scaffold length,min S scaffold length,max S kmer coverage,min S kmer coverage,avg S scaffold pairwise distance,min S scaffold pairwise distance,min S scaffold pairwise distance" > S_all_scaffold_stats.csv

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
  echo "${sample},s,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> S_all_scaffold_stats.csv;
done



#filter M segment matches!
#print output to file with the following headers: 
echo "sample,segment,num M scaffolds,max M scaffold length,min M scaffold length,max M kmer coverage,min M kmer coverage,avg M scaffold pairwise distance,min M scaffold pairwise distance,min M scaffold pairwise distance" > M_all_scaffold_stats.csv

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
  echo "${sample},m,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> M_all_scaffold_stats.csv;
done



#L segment matches filtering!
#print output to file with the following headers: 
echo "sample,segment,num L scaffolds,max L scaffold length,min L scaffold length,max L kmer coverage,min L kmer coverage,avg L scaffold pairwise distance,min L scaffold pairwise distance,min L scaffold pairwise distance" > L_all_scaffold_stats.csv

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
  echo "${sample},l,$(grep -c '>' ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds.fasta),${sort_filt_len_values[0]},${sort_filt_len_values[-1]},${sort_filt_cov_values[0]},${sort_filt_cov_values[-1]}" >> L_all_scaffold_stats.csv;
done



#compute pairwise identity with clustalw
#sed -i '1{s/$/,avg scaffold pairwise identity,/}' L_all_scaffold_stats.csv


for file in *${sample_suffix}_output/*${sample_suffix}_L_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
   
  cp ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds.fasta test;
  clustalw test;clustalw -outputtree=dist -tree -infile=test.aln;l test.aln;lls test.dst;
  cp test.dst ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_clustal.dst;
  dist_arr=();
  sort_dist_arr=();
  sum=0;
  for line in $(grep 'NODE' test.dst); 
    do for word in $line; 
      do if (( $(echo "$word > 0" |bc -l) )); then
        echo $word;
        dist_arr+=("${word}");
        sum=$(echo "$sum + $word"|bc -l);
      fi;
    done;
  done
len=${#dist_arr[@]}
sort_dist_arr=(`printf '%s\n' "${dist_arr[@]}" | sort -gr`); 
avg=$(echo "$sum/$len"|bc -l)
line=$(grep -no "${sample}," L_all_scaffold_stats.csv);

sed -i "${line%%:*}{s/$/,${avg},${sort_dist_arr[0]},${sort_dist_arr[-1]}/}" L_all_scaffold_stats.csv;

rm test*
done



for file in *${sample_suffix}_output/*${sample_suffix}_M_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
   
  cp ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds.fasta test;
  clustalw test;clustalw -outputtree=dist -tree -infile=test.aln;l test.aln;lls test.dst;
  cp test.dst ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_clustal.dst;
  dist_arr=();
  sort_dist_arr=();
  sum=0;
  for line in $(grep 'NODE' test.dst); 
    do for word in $line; 
      do if (( $(echo "$word > 0" |bc -l) )); then
        echo $word;
        dist_arr+=("${word}");
        sum=$(echo "$sum + $word"|bc -l);
      fi;
    done;
  done
len=${#dist_arr[@]}
sort_dist_arr=(`printf '%s\n' "${dist_arr[@]}" | sort -gr`); 
avg=$(echo "$sum/$len"|bc -l)
line=$(grep -no "${sample}," M_all_scaffold_stats.csv);

sed -i "${line%%:*}{s/$/,${avg},${sort_dist_arr[0]},${sort_dist_arr[-1]}/}" M_all_scaffold_stats.csv;

rm test*
done

for file in *${sample_suffix}_output/*${sample_suffix}_S_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
   
  cp ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds.fasta test;
  clustalw test;clustalw -outputtree=dist -tree -infile=test.aln;l test.aln;l test.dst;
  cp test.dst ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_clustal.dst;
  dist_arr=();
  sort_dist_arr=();
  sum=0;
  for line in $(grep 'NODE' test.dst); 
    do for word in $line; 
      do if (( $(echo "$word > 0" |bc -l) )); then
        echo $word;
        dist_arr+=("${word}");
        sum=$(echo "$sum + $word"|bc -l);
      fi;
    done;
  done
len=${#dist_arr[@]}
sort_dist_arr=(`printf '%s\n' "${dist_arr[@]}" | sort -gr`); 
avg=$(echo "$sum/$len"|bc -l)
line=$(grep -no "${sample}," S_all_scaffold_stats.csv);

sed -i "${line%%:*}{s/$/,${avg},${sort_dist_arr[0]},${sort_dist_arr[-1]}/}" S_all_scaffold_stats.csv;

rm test*
done




#get reads aligned to L segment assembly with bowtie2: 

for file in *${sample_suffix}_output/*${sample_suffix}_L_covlenfilt_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
  bowtie2-build ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds.fasta ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds_DB;
  bowtie2 -p 8 -x ${sample}${sample_suffix}_output/${sample}${sample_suffix}_L_covlenfilt_scaffolds_DB -1 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R1.qhtrim.unmapped.sorted2.fastq  -2 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R2.qhtrim.unmapped.sorted2.fastq  -S ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_L_covlenfilt_assembly.sam > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_L_covlenfilt_assembly_bowtie2log.txt 2>&1;
  match=$(grep "aligned concordantly exactly 1 time" ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_L_covlenfilt_assembly_bowtie2log.txt);
  matchtemp=${match%%%*};
  echo;
  echo "match for aligned concordantly output line = ${match}, and after 1st parse step matchtemp = ${matchtemp}";
  concordant=${matchtemp##*\(};
  echo;
  echo "% concordant reads variable currently set to ${concordant} for sample ${sample}";
  linematch=$(grep -no "${sample}," L_all_scaffold_stats.csv);
  sed -i "${linematch%%:*}{s/$/,${concordant}/}" L_all_scaffold_stats.csv;
  done




#get reads aligned to M segment assembly with bowtie2: 

for file in *${sample_suffix}_output/*${sample_suffix}_M_covlenfilt_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
  bowtie2-build ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds.fasta ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds_DB;
  bowtie2 -p 8 -x ${sample}${sample_suffix}_output/${sample}${sample_suffix}_M_covlenfilt_scaffolds_DB -1 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R1.qhtrim.unmapped.sorted2.fastq  -2 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R2.qhtrim.unmapped.sorted2.fastq  -S ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_M_covlenfilt_assembly.sam > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_M_covlenfilt_assembly_bowtie2log.txt 2>&1;
  match=$(grep "aligned concordantly exactly 1 time" ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_M_covlenfilt_assembly_bowtie2log.txt);
  matchtemp=${match%%%*};
  echo;
  echo "match for aligned concordantly output line = ${match}, and after 1st parse step matchtemp = ${matchtemp}";
  concordant=${matchtemp##*\(};
  echo;
  echo "% concordant reads variable currently set to ${concordant} for sample ${sample}";
  linematch=$(grep -no "${sample}," M_all_scaffold_stats.csv);
  sed -i "${linematch%%:*}{s/$/,${concordant}/}" M_all_scaffold_stats.csv;
  done





#get reads aligned to S segment assembly with bowtie2: 

for file in *${sample_suffix}_output/*${sample_suffix}_S_covlenfilt_scaffolds.fasta; 
  #get filename
  do fileonly=$(basename $file); 
  #get sample letter from file namel
  #sample=${fileonly##*/};
  sample=${fileonly%%_*};
  bowtie2-build ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds.fasta ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds_DB;
  bowtie2 -p 8 -x ${sample}${sample_suffix}_output/${sample}${sample_suffix}_S_covlenfilt_scaffolds_DB -1 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R1.qhtrim.unmapped.sorted2.fastq  -2 /group/sldmunozgrp/cm_cysto_miseq_M1382P_Mattson/final_trimmed_fastqs/${sample}_SXX_L001_001.R2.qhtrim.unmapped.sorted2.fastq  -S ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_S_covlenfilt_assembly.sam > ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_S_covlenfilt_assembly_bowtie2log.txt 2>&1;
  match=$(grep "aligned concordantly exactly 1 time" ${sample}${sample_suffix}_output/${sample}${sample_suffix}_reads_to_S_covlenfilt_assembly_bowtie2log.txt);
  matchtemp=${match%%%*};
  echo;
  echo "match for aligned concordantly output line = ${match}, and after 1st parse step matchtemp = ${matchtemp}";
  concordant=${matchtemp##*\(};
  echo;
  echo "% concordant reads variable currently set to ${concordant} for sample ${sample}";
  linematch=$(grep -no "${sample}," S_all_scaffold_stats.csv);
  sed -i "${linematch%%:*}{s/$/,${concordant}/}" S_all_scaffold_stats.csv;
  done

  
