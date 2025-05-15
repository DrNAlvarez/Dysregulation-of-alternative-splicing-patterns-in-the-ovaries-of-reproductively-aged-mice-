# Dysregulation of alternative splicing patterns in the ovaries of reproductively aged mice 

Aging is associated with widespread changes in cellular splicing patterns, which alter gene isoform expression and function. Splicing changes are increasingly recognized as contributors to aging phenotypes, rather than mere byproducts of declining RNA quality control. Ovarian aging, a key aspect of female reproductive decline, involves a reduction in both follicle number and function, leading to diminished estrogen production and the onset of menopause. This process also compromises oocyte quality, elevating the risks of genetic disorders, miscarriage, and infertility. Despite its significance, how splicing landscapes shift within the aging ovary remains poorly understood. To investigate the role of alternative splicing in ovarian aging, we performed direct RNA sequencing (dRNA-seq) on ovaries from reproductively young and aged mice, capturing the transcriptome at single-molecule, full-length resolution. Comparative analysis revealed that reproductive aging is marked by increased intron retention and significant changes in exon usage, particularly within the first and last exons. Intron retention, a common feature in aging tissues, is often linked to impaired protein expression, suggesting potential disruptions in gene function. Additionally, we observed a shift toward more distal alternative last exons and more proximal alternative first exons in aged ovaries, leading to the production of longer transcripts. These findings imply that aging may alter promoter utilization and transcription termination, possibly through shifts in polyadenylation site selection. Our work underscores the potential of alternative splicing as a driver of ovarian aging and reveals new avenues for understanding the molecular underpinnings of reproductive decline. By leveraging the power of direct RNA sequencing, we provide insights into novel protein isoforms that may play crucial roles in ovarian physiology and reproductive aging.

# Usage
Steps

1. Aggregate all fastq data and run flair (https://github.com/BrooksLabUCSC/flair) to generate a single gtf_reference file
   
    zcat *.fastq.gz > all.fastq.gz

    flair align -g GRCm38_68.fa -r all.fastq.gz  -o all.flair.aligned -t 120 --nvrna

    flair correct --nvrna -g GRCm38_68.fa -q all.flair.aligned.bed -f Mus_musculus.GRCm38.102.gtf -o all.correct.flair -t 120

    flair collapse -g GRCm38_68.fa -r all.fastq.gz -q all.correct.flair_all_corrected.bed -f Mus_musculus.GRCm38.102.gtf -o all.flair.collapse --generate_map --annotation_reliant generate -s 1 --temp_dir flair_temp/ --keep_intermediate --check_splice --stringent -t 120

    flair quantify -r reads_manifest_all.tsv -i all.flair.collapse.isoforms.fa --generate_map --stringent --check_splice --isoform_bed all.flair.collapse.isoforms.bed --tpm -t 120 --temp_dir tmp/ -o coc_all

2. Run bioseq2seq to translate cDNA into proteins (https://github.com/josephvalencia/bioseq2seq). The output will need to be parsed into a fasta file
   python ~/bioseq2seq/bin/translate.py --checkpoint ~/bioseq2seq/best_bioseq2seq_mammalian_200-1200.pt --input all.flair.collapse.isoforms.fa --o all.flair.collapse.translation --mode bioseq2seq --num_gpus 1 --beam_size 1 --n_best 1 --max_tokens 1200 --max_decode_len 400
   
   python parse_bioseq2seq.py -i all.flair.collapse.isoforms_full_preds.txt -o parsed_sequences


For Figure 1.

Run the splicing annotation script to make a .gtf file with the splicing events annnotated. This will generate files needed for visualization and analysis:
   
    all.flair.collapse.isoforms_event_label_promoter_group_data.pkl
    all.flair.collapse.isoforms_event_label.gtf
    ASevents.csv

    python annotate_gtf_with_splicing.py -o all.flair.collapse.isoforms -p 120 Mus_musculus.GRCm38.102.gtf all.samples.tpm.tsv


Run fetch_ensembl_biotype.py to get biotypes for the transcripts annotated in all.flair.collapse.isoforms_event_label.gtf

    python fetch_ensembl_biotype.py --host http://useastdb.ensembl.org --dataset mmusculus_gene_ensembl -o ensembl_mouse_ids.tsv

Run translate_cDNA.py to translate all cDNA identified. 

    python check_cdna_translation.py all.flair.collapse.isoforms.pep/parsed_sequences.fasta all.flair.collapse.isoforms.fa translation_start_stops.tsv -p 120

You will now be able to run Figure_1 notebook

For Figure 2. 

Run run_diff_expr.py  to generate a data frame of differentially spliced isoforms between two conditions only. If you have more conditions to compare you will have run them in pairs. (A vs B or A vs C ect. For the manuscript we have two condtions a_vs_c (YNG PMSG vs AMA PMSG) and b_vs_d (YNG PMSG+hCG vs AMA PMSG+hCG)

    run_diff_expr.py --asevents ASevents_ovary_v4.csv --counts all.samples.counts.tsv --condA 1 2 3 --condB 7 8 9 --groups gene_id promoter_group --nsamples 200 --ntune 200 --cores 8 --out results_splice.csv

You will now be able to run Figure_2 notebook, it requires the following files

    all.flair.collapse.isoforms.gtf
    results_splice.csv (a_vs_c and b_vs_d)
    ASevents.csv
    
For Figure 3. 

You will need the following files:

    all.flair.collapse.isoforms_event_label.gtf
    results_splice.csv (a_vs_c and b_vs_d)
    ASevents.csv

For Figure 4.

You will need to run the following scripts to generate the necessary files:

    python translate_cDNA.py  parsed_sequences.fasta all.flair.collapse.isoforms.fa translation_start_stops.csv -p 96
    
    python add_orf_features.py --gtf all.flair.collapse.isoforms_event_label.gtf --codon translation_start_stops.csv --out all.flair.collapse.isoforms_event_label_start_stop_orf.gtf --procs 96

You will need the following files as input for Figure 4:

    results_splice.csv (a_vs_c and b_vs_d)
    ASevents.csv
    all.flair.collapse.isoforms_event_label_start_stop_orf.gtf
    all.flair.collapse.isoforms_event_label_promoter_group_data.pkl

For Figure 5. 


You will need to run the following code to generate the necessary files
   
    python build_translation_promoter_files.py --translation translation_start_stops.csv --gtf all.flair.collapse.isoforms_event_label_start_stop_orf.gtf --csv translation_data_by_promoter_group.csv --fasta translation_data_by_promoter_group.fa

    hmmpress Pfam-A.hmm

    hmmscan --domtblout output.domtblout Pfam-A.hmm translation_data_by_promoter_group.fa

    python annotate_pfam_domains.py --domtbl output.domtblout --translation translation_start_stops.csv --gtf all.flair.collapse.isoforms_event_label.gtf --out-prefix isoforms_event_label --evalue 1e-5 --procs 96

    output files are:
    isoforms_event_label.pfam_start_stops.tsv
    isoforms_event_label.pfam_start_stop_with_translation.tsv
    isoforms_event_label_start_stop_pfam.gtf

You will need the following files as input for Figure 5

    all.flair.collapse.isoforms_event_label_start_stop_orf.gtf
    all.flair.collapse.isoforms_event_label_start_stop_pfam.gtf (isoforms_event_label_start_stop_pfam.gtf)
    all.flair.collapse.isoforms_event_label.gtf
    pfam_start_stops.tsv (isoforms_event_label.pfam_start_stops.tsv)
    results_splice.csv (a_vs_c and b_vs_d)
    ASevents.csv
    shortNdufs4.png (alpha fold image output)
    Full_lengthNdufs4.png (alpha fold image output)
