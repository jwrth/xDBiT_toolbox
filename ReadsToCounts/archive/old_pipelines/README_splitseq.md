# Split-seq toolbox

 This is a toolbox for processing raw sequencing output from Split-seq experiments into a digital gene expression matrix that will contain integer counts of the number of transcripts per single cell. It provides a bash script that serves as a wrapper for multiple analysis steps, including demultiplexing of the raw data by molecular (UMI) and cellular barcode, filtering cellular barcodes by a list of expected barcodes, alignment of reads to a reference genome, collecting basic QC metrics and counting UMIs per cell.

 The pipeline uses some custom python scripts, and many tools from the [Drop-seq toolbox](https://github.com/broadinstitute/Drop-seq/releases) (Mc Caroll lab, Harvard Medical school) as well as [Picard](https://broadinstitute.github.io/picard/) (Broad institute), which are all included in this toolbox.

 Please refer to the [user manual](./man/Splitseq_toolbox_manual.pdf) for a description of how to use it.

## Contact
This toolbox is provided by the [Snijder lab](https://www.snijderlab.org) at ETH Zurich.

If you have questions or find a bug, son't hesitate to contact Rebekka Wegmann: [wegmann@imsb.biol.ethz.ch](mailto:wegmanre@imsb.biol.ethz.ch)

## Notes and caution
This tool comes with no warranty and accurate function is not guaranteed. It laso does not exactly recapitulate the approach described in the SPLiT-seq paper.
Currently, we use only the poly-T RT primers, therefore **this tool does not collapse the random hexamer and polyT barcodes!**
