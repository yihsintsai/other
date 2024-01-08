#!/bin/bash


for f in $(ls -v ./|grep ONT) ; do
        sample_id=$( echo ${f} | sed -e 's/^.*A/A/g' |sed -e 's/.tar$//g' |uniq)
                id=$( echo ${f}|sed -e 's/.tar$//g' |uniq)
        if [[ $(echo ${f} | tr "_" "\t" | cut -f 2 |sed -e 's/.tar$//g') != "${sample_id}" ]]; then
                batch=$(ls ${f} | tr "_" "\t" | cut -f 2)
                sample=$(echo ${sample_id}"_"${batch})
                echo $batch
                echo "sample =" ${sample}
        else
                echo "sample =" ${sample}
        fi

        sshpass -p "10216116irislove" \
        rsync -avP ./${f} -e \
        ssh yh@172.27.149.102:/mnt/nas2/GBM/ONT/${sample}/
        wait
        echo $sample_id "upload done" >> ${id}.txt
        sshpass -p $(cat pass.txt) \
        rsync -avP ./${id}.txt -e \
        ssh yh@172.27.149.102:/mnt/nas2/GBM/ONT/${sample}/
done &
