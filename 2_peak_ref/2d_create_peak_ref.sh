#!/usr/bin/bash

## This script creates peak ref using MACS2 callpeak outputs on polyA reads and such polyA reads ##

# Author: Surbhi Sona, Ninh Le

## User inputs ##

OPTIND=1

while getopts "m:p:o:s:c:n:t:" opt;do
  case $opt in
  m)
    macs_dir=$OPTARG
    ;; #-m path to MACS2 peaks output
  p)
    polya_reads_dir=$OPTARG
    ;; #-p path to polya dir
  o)
    out_dir=$OPTARG
    ;; #-o specifies output directory for peak ref
  s) 
    strand=$OPTARG
    ;; #-s specify strand as plus or minus
  c)
    chr=$OPTARG
    ;; # -c specify chr (1-22 or X or Y or multiple chromosomes separated by _)
  n)
    min_polya=$OPTARG
    ;; #-n minimum polya reads that a peak must be supported by
   \?) echo "Option ${opt} not valid";;   
  t)
    script_dir=$OPTARG
    ;; #-t script dir path 
  esac
done

function execute_command {
    echo "$(date +"%Y-%m-%d %H:%M:%S"):"
    "$@"
}

### Set paths ###

cd ${out_dir}
outdir=${out_dir}${strand}/
mkdir ${outdir}

peaks_per_chr=${outdir}0_peaks_per_chr/
polya_per_chr=${outdir}0_polya_per_chr/
int_dir=${outdir}1_int_dir/
peaks_count_dir=${outdir}2_peaks_count_dir/
peak_ref_table=${outdir}3a_peak_ref_table/
pr_dir=${outdir}3b_pr_dir/
final_peak_ref=${outdir}4_polya_supported_peak_ref/

### 0. Prepare files ###

#(a) create dir

cd ${outdir}

mkdir -p ${peaks_per_chr}
mkdir -p ${polya_per_chr}
mkdir -p ${int_dir}
mkdir -p ${peaks_count_dir}
mkdir -p ${peak_ref_table}
mkdir -p ${pr_dir}
mkdir -p ${final_peak_ref}
 

# Create chr file (currently doing it manually to exclude non standard chr)

chr=$(echo $chr | tr '_' ' ')

# (b) Specify peaks and polya  files

cd ${macs_dir}
peaks_file=${macs_dir}*${strand}*peaks.narrowPeak
polya_file=${polya_reads_dir}*${strand}*nowrongstrand.bed

# (c) Split peak columns by chr, subset desired col and add strand column

execute_command echo Creating peak from $peaks_file and polya files from $polya_file by chr: ${chr}

for c in $(echo $chr); do
execute_command echo Creating peak ref for ${c}...

if [ ${strand} == 'plus' ]
then
   cat ${peaks_file} | grep -w ${c} | awk 'OFS="\t" {print $1,$2,$3,$4,$5="+"}' > ${peaks_per_chr}peaks_${c}_${strand}.bed
else
   cat ${peaks_file} | grep -w ${c} | awk 'OFS="\t" {print $1,$2,$3,$4,$5="-"}' > ${peaks_per_chr}peaks_${c}_${strand}.bed
fi
 
# subset poly table too

cat ${polya_file} | grep -w ${c} | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > ${polya_per_chr}polya_${c}_${strand}.bed


# (e) sort 

execute_command echo Sort peaks by start and end coordinates... 
cat  ${peaks_per_chr}peaks_${c}_${strand}.bed | sort -k2,2n -k3,3n >  ${peaks_per_chr}peaks_${c}_sorted_${strand}.bed
execute_command echo Sort polyA reads by start and end coordinates...
cat  ${polya_per_chr}polya_${c}_${strand}.bed | sort -k2,2n -k3,3n >  ${polya_per_chr}polya_${c}_sorted_${strand}.bed

### 1a. Run intersect ###

execute_command echo Intersect...

bedtools intersect -a ${peaks_per_chr}peaks_${c}_sorted_${strand}.bed -b ${polya_per_chr}polya_${c}_sorted_${strand}.bed -wa -wb -loj > ${int_dir}peaks_int_polya_${c}_${strand}.bed

### 1b. Sort intersect table ###

execute_command echo Sort intersect table...
cat ${int_dir}peaks_int_polya_${c}_${strand}.bed | sort -k6,6  > ${int_dir}peaks_int_polya_${c}_sorted_${strand}.bed

### 1c. Tag non-poly peaks

execute_command echo Tag polyA peaks...
awk  'FS=OFS="\t" {if($6==".") {$12="non_polya"} else {$12="polya"}; {print $0}}' ${int_dir}peaks_int_polya_${c}_sorted_${strand}.bed > ${int_dir}peaks_int_polya_${c}_tagged_${strand}.bed

### 1d. Sort by polya column
cat ${int_dir}peaks_int_polya_${c}_tagged_${strand}.bed | sort -k12,12 > ${int_dir}peaks_int_polya_${c}_tagged_sorted_${strand}.bed

### 2. Create peaks count

execute_command echo Count how many polyA reads each polyA intersect with and sort peaks by count...
cat ${int_dir}peaks_int_polya_${c}_tagged_${strand}.bed | grep -w polya | awk 'FS=OFS="\t" {print $4}' | sort | uniq -c | sort -k1,1n > ${peaks_count_dir}peaks_${c}_count_sorted_${strand}.txt

### 3a. Create peak ref table ###
execute_command echo Create peak ref table from ${int_dir}peaks_int_polya_${c}_tagged_sorted_${strand}.bed...
cat ${int_dir}peaks_int_polya_${c}_tagged_sorted_${strand}.bed | sort -u -k4,4  > ${peak_ref_table}peak_ref_${c}_${strand}.bed

### 3b. Create inferred PR table ### (contains 3 columns - peakID, inf_pr_start and inf_pr_end)

f=$(awk '{print $2}' ${peaks_count_dir}peaks_${c}_count_sorted_${strand}.txt)
execute_command echo Create peak ref table with the inferred PR, which are the ranges that encompass the polyA sites of all the intersected polyA reads...
if test -f ${pr_dir}pr_${c}_${strand}.bed; then echo Delete the existing ${pr_dir}pr_${c}_${strand}.bed; rm ${pr_dir}pr_${c}_${strand}.bed; fi

for peaks in ${f[*]};do
cat ${int_dir}peaks_int_polya_${c}_sorted_${strand}.bed | grep -w ${peaks}  | awk '{if(min==""){min=max=$7}; if($8>max) {max=$8}; if($7< min) {min=$7};} END {print $4"\t"min"\t"max}' >> ${pr_dir}pr_${c}_${strand}.bed
done

# 3c. Sort
execute_command echo Sort peaks by the starts and ends of the inferred PRs...
cat ${pr_dir}pr_${c}_${strand}.bed | sort -k2,2n -k3,3n > ${pr_dir}pr_${c}_sorted_${strand}.bed

### 3d. Add inferred PR, peak & PR width and polyA count  ###
execute_command echo Process the last peak ref table, including filtering for peaks with at least ${min_polya} reads
Rscript ${script_dir}2d_update_peak_ref.R ${outdir} ${strand} ${c} ${min_polya}

execute_command echo Finish creating peak ref for ${c}!

done
