# Dysregulation of alternative splicing patterns in the ovaries of reproductively aged mice 

Aging is associated with widespread changes in cellular splicing patterns, which alter gene isoform expression and function. Splicing changes are increasingly recognized as contributors to aging phenotypes, rather than mere byproducts of declining RNA quality control. Ovarian aging, a key aspect of female reproductive decline, involves a reduction in both follicle number and function, leading to diminished estrogen production and the onset of menopause. This process also compromises oocyte quality, elevating the risks of genetic disorders, miscarriage, and infertility. Despite its significance, how splicing landscapes shift within the aging ovary remains poorly understood. To investigate the role of alternative splicing in ovarian aging, we performed direct RNA sequencing (dRNA-seq) on ovaries from reproductively young and aged mice, capturing the transcriptome at single-molecule, full-length resolution. Comparative analysis revealed that reproductive aging is marked by increased intron retention and significant changes in exon usage, particularly within the first and last exons. Intron retention, a common feature in aging tissues, is often linked to impaired protein expression, suggesting potential disruptions in gene function. Additionally, we observed a shift toward more distal alternative last exons and more proximal alternative first exons in aged ovaries, leading to the production of longer transcripts. These findings imply that aging may alter promoter utilization and transcription termination, possibly through shifts in polyadenylation site selection. Our work underscores the potential of alternative splicing as a driver of ovarian aging and reveals new avenues for understanding the molecular underpinnings of reproductive decline. By leveraging the power of direct RNA sequencing, we provide insights into novel protein isoforms that may play crucial roles in ovarian physiology and reproductive aging.

# Usage

A GTF file generated from FLAIR pipeline analysis is required. https://github.com/BrooksLabUCSC/flair
Below is the code used for our use case, your requirements may vary.

Note: Directories have been omitted for security reseaons and replaced with generic names.

Note: We joined all the fastq files to generate a experiment specific GTF. 

zcat path/to/fastqs/*.fastq.gz  > all.fastq.gz

flair align -g GRCm38_68.fa -r all.fastq.gz  -o all.flair.aligned -t 120 --nvrna

flair correct --nvrna -gGRCm38_68.fa -q all.flair.aligned.bed -f Mus_musculus.GRCm38.102.gtf -o all.flair -t 120

flair collapse -g GRCm38_68.fa -r all.fastq.gz -q all.correct.flair_all_corrected.bed -f Mus_musculus.GRCm38.102.gtf -o all.flair.collapse --generate_map --annotation_reliant generate -s 1 --temp_dir flair_temp/ --keep_intermediate --check_splice --stringent -t 120

flair quantify -r reads_manifest_all.tsv -i all.flair.collapse.isoforms.fa --generate_map --stringent --check_splice --isoform_bed all.flair.collapse.isoforms.bed --tpm -t 48 --temp_dir tmp/ -o all.flair
