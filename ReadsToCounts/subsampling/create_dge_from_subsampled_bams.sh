#!/usr/bin/env bash

## Create DGE matrices from subsampled files

dropseq_root="../external_tools/Drop-seq_tools-2.1.0/"
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

outdir=`pwd`
tmpdir="tmp"
mkdir $tmpdir

files=gene_function_tagged_*

for file in $files
do
	echo "Processing of ${file} started..."
	filename="$(basename ${file} .bam)"
	perc=`echo ${filename} | awk -F "_" '{print $NF}'`
	
	# counting exonic reads only
	dge="${dropseq_root}/DigitalExpression I=${outdir}/${file} O=${outdir}/DGE_matrix_${perc}.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 TMP_dir=${tmpdir}"
	$echo_prefix $dge

	# counting both intronic and exonic reads
	dge_with_introns="${dropseq_root}/DigitalExpression I=${outdir}/${file} O=${outdir}/DGE_matrix_with_introns_${perc}p.txt.gz READ_MQ=10 EDIT_DISTANCE=1 MIN_NUM_GENES_PER_CELL=100 LOCUS_FUNCTION_LIST=INTRONIC TMP_dir=${tmpdir}"
	$echo_prefix $dge_with_introns

	# collect RNAseq metrics with PICARD
	#rnaseq_metrics="java -jar ${picard_jar} CollectRnaSeqMetrics I=${outdir}/${file} O=${outdir}/rnaseq_metrics_${perc}.RNA_Metrics REF_FLAT=${refflat} STRAND=FIRST_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=${rRNA_intervals}"
	#$echo_prefix $rnaseq_metrics

	echo "Processing of ${file} finished."
done

echo "Processing of all files finished."