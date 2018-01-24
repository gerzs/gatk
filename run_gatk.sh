for i in {01..48}
do
        cp /tusinas/databases/tcga/alignment/LIHC/LIHC_${i}_01/star.Aligned.sortedByCoord.out.bam /var/tmp/snp/${i}_tumor.bam &
        pid1=$!
        cp /tusinas/databases/tcga/alignment/LIHC/LIHC_${i}_11/star.Aligned.sortedByCoord.out.bam /var/tmp/snp/${i}_normal.bam &
        pid2=$!

        for j in tumor normal
        do
                if [[ $j == 'tumor' ]]
                then
                        wait $pid1
                        samtools view -H ${i}_tumor.bam > ${i}_tumor.header.sam
                        sed 's/SO:sorted//' ${i}_tumor.header.sam  | sed "s/\t\t/\t/" > ${i}_tumor.header_corrected.sam
                        samtools reheader ${i}_tumor.header_corrected.sam ${i}_tumor.bam > ${i}_tumor.newheader.bam
                        java -jar picard.jar AddOrReplaceReadGroups I=${i}_tumor.newheader.bam O=${i}_tumor.rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
                        java -jar picard.jar MarkDuplicates I=${i}_${j}.rg_added_sorted.bam O=${i}_${j}.dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${i}_${j}.output.metrics
                        java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R /bigdisk/databases/genomes/human/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -I ${i}_${j}.dedupped.bam -o ${i}_${j}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALL &
                        pid3=$!
                else
                        wait $pid2
                        samtools view -H ${i}_${j}.bam > ${i}_${j}.header.sam
                        sed 's/SO:sorted//' ${i}_${j}.header.sam  | sed "s/\t\t/\t/" > ${i}_${j}.header_corrected.sam
                        samtools reheader ${i}_${j}.header_corrected.sam ${i}_${j}.bam > ${i}_${j}.newheader.bam
                        java -jar picard.jar AddOrReplaceReadGroups I=${i}_${j}.newheader.bam O=${i}_${j}.rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
                        java -jar picard.jar MarkDuplicates I=${i}_${j}.rg_added_sorted.bam O=${i}_${j}.dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${i}_${j}.output.metrics
                        java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R /bigdisk/databases/genomes/human/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -I ${i}_${j}.dedupped.bam -o ${i}_${j}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALL &
                        pid4=$!
                fi
        done
        wait $pid3 $pid4
        java -jar GenomeAnalysisTK.jar -T MuTect2 -R /bigdisk/databases/genomes/human/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -I:tumor ${i}_tumor.split.bam -I:normal ${i}_normal.split.bam --dbsnp /bigdisk/users/gerzs/snptest/All_20170710.vcf --cosmic /bigdisk/users/gerzs/snptest/CosmicCodingMuts.vcf -o /bigdisk/users/gerzs/snptest/LIHC/LIHC_${i}.vcf -nct 2 -U ALL &

        if [[ -f /bigdisk/users/gerzs/snptest/LIHC/LIHC_${i}.vcf.idx ]]
        then
                rm ${i}_tumor.dedupped.bai
                rm ${i}_normal.dedupped.bai
                rm ${i}_tumor.dedupped.bam
                rm ${i}_normal.dedupped.bam
                rm ${i}_tumor.header_corrected.sam
                rm ${i}_normal.header_corrected.sam
                rm ${i}_tumor.header.sam
                rm ${i}_normal.header.sam
                rm ${i}_tumor.newheader.bam
                rm ${i}_normal.newheader.bam
                rm ${i}_tumor.output.metrics
                rm ${i}_normal.output.metrics
                rm ${i}_tumor.rg_added_sorted.bam
                rm ${i}_normal.rg_added_sorted.bam
                rm ${i}_tumor.split.bam
                rm ${i}_tumor.split.bai
                rm ${i}_normal.split.bam
                rm ${i}_normal.split.bai
        fi
done
