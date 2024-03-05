#!/bin/bash

for i in $(cat 55315445_region) 
#for i in 00162a68-c663-453e-b57a-7db279c3e12a
	do
    echo ${i} >> read_info
		grep ${i} ./A23.bed > seq.read
	  samtools view -@ 12 /mnt/d/GBM/3.ONT/A23.fq.gz.candidates.bam | grep ${i} |cut -f 3,4,6,10 >  tem_seq.txt
		grep ${i} ../A23.fq.gz.candidates.fastq -A3 | sed -n 2p > fa_seq_tmp.txt

    ##check str number 
    str_num=$(grep ${i} region_s_tail.txt | cut -f 2 )
    echo "str_position" ${str_num} >> read_info
    read_num=$(grep -n ${str_num} seq.read | awk -F: '{print $1}')
    echo "which line" ${read_num} >> read_info
    if [[ $(echo ${read_num} |wc -w ) > 1  ]] ; then
      echo -e ${i}'\t'${read_num}'\t'"multiple maping, need check " >> q_reads
    else
    total_num=$(cut -f 2 seq.read  |wc -l)
    glo_num=$((${read_num}+1))
    if [[ $total_num < ${glo_num} ]] ; then
      echo no H clip" >> read_info
    else    
      strand=$(sed -n "${glo_num}p" seq.read| cut -f 6 )
      echo "F or R strand" ${strand} >> read_info
      read_position=$(sed -n "${glo_num}p" seq.read| cut -f 2 )
      final_num=$((${read_position}+1))
      echo "global read str position" ${final_num} >> read_info
      wait
      if [[ ${strand} == "+" ]] ; then
        echo "For strand" >> read_info
        seq=$(grep ${final_num} tem_seq.txt | cut -f 4 )
        seq_count=$(echo ${seq}|wc -c)
        echo "total amount" ${seq_count} >> read_info
        if [[ (${seq_count} > 2000) ]] ; then
            seq_part=$(echo ${seq} | cut -c 1-2000)
        else 
            seq_part=${seq}
            fi
        H_num=$(sed -n "${glo_num}p" seq.read | cut -f 7|grep -oe "^[0-9]*H" |sed -e 's/H//g')
        echo "H-clip number" ${H_num} >> read_info
        seq_position=$(grep -bo ${seq_part} fa_seq_tmp.txt | awk -F: '{print $1}')
        echo "read str" ${seq_position} >> read_info
        H_position=$((${seq_position}-${H_num}))
        if [[ ${H_position} == 0  ]] ; then
          H_position_str=$((${H_position}+1))
          echo "H-clip str" ${H_position} +1  >> read_info
        elif [[ ${H_position} < 0 ]] ; then
          H_position_str=${H_position}
        else
          echo ${i} "sequence err" >>  read_info
        fi  
          H_seq=$(cut -c ${H_position_str}-${seq_position} fa_seq_tmp.txt)
          echo -e ${i}"_H"'\t'${H_num}"H"'\t'${H_seq} >>  seq_20240224.txt
      else 
        echo "reverse strand" >> read_info
        seq=$(echo $(grep ${final_num} tem_seq.txt | cut -f 4 ) | tr 'ATCGatcg' 'TAGCtagc' | rev)
        seq_count=$(echo ${seq}|wc -c)
        echo "total amount" ${seq_count} >> read_info
        if [[ (${seq_count} > 2000) ]] ;then
            seq_part=$(echo ${seq} | cut -c 1-2000)
        else 
            seq_part=${seq}
            fi
        H_num=$(sed -n "${glo_num}p" seq.read | cut -f 7 | grep -oe "^[0-9]*H" |sed -e 's/H//g')
        echo "H-clip number" ${H_num} >> read_info
        seq_position=$(grep -bo ${seq_part} fa_seq_tmp.txt | awk -F: '{print $1}')
        echo "read str position" ${seq_position} >> read_info
        seq_end_pos=$((${seq_position}+${seq_count}))
        echo "read end position" ${seq_end_pos} >> read_info
        H_position=$((${seq_end_pos} + ${H_num}))
        echo "read end position" ${H_position} >> read_info
        total_fa_num=$(wc -c fa_seq_tmp.txt |cut -f 1 -d ' ' )
        echo 
        if [[ ${total_fa_num} <  ${H_position} ]] ; then
          echo ${i} "sequence err" >> read_info
        else
          H_position_end=${H_position}
          fi
        #H_seq=$(cut -c ${seq_end_pos}-${H_position_end} fa_seq_tmp.txt | tr 'ATCGatcg' 'TAGCtagc' | rev )
        H_seq=$(cut -c ${seq_end_pos}-${H_position_end} fa_seq_tmp.txt  )
        echo -e ${i}"_H"'\t'${H_num}"H"'\t'${H_seq} >>  seq_20240226.txt
        fi
      fi
    fi 
    wait
        
      done
