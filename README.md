These are the scripts used for the manuscript: **Mapping DNA Methylation to Cardiac Pathologies Induced by Beta-Adrenergic Stimulation in a Large Panel of Mice**  Doi: https://doi.org/10.1080/15592294.2025.2524411

Script 1:  Methylation_Calling_Script.sh - BASH

  Input:  fastQ.gz files from the sequencer, fasta files for the whole genome sequence
  Function: It constructs an RRBS_digested genome, followed by BSseeker2 alignment, followed by CG calling and SNP identification
  Output:  One file per input with CGs and one file per input of identified SNPs

  Script 2: ProcessCGmap_COMBAT_Filter.R - R
    Input:  Outputs from Script 1, Sample Info with Batch Effects, Choice of coverage filters
    Function:  Combines all files into one, runs COMBAT to eliminate batch effects, filters based on chosen filters (coverage/ % not NA), and calculates hypervariable CpGs
    Output:  A number of files containing the merged CpG methylation statuses across the entire cohort at various stages (pre COMBAT, no filters, only hypervariable, etc)

  Script 3: RunMacau.R - R
    Input: Hypervariable Methylation and Coverage files from step 2, Phenotype Data, 
    Function:  Separates data into control/ISO subsets, calculates kinship matrices, and runs MACAU.  Combines all MACAU runs together, generates Manhattan and QQ plots
    Output:  Individual files and merged pvalue file, manhattan/QQ plots for every phenotype

  Script 4: LocusZoomMouse.R - R
    Function:  The Locuszoomr package only works with human data.  This function allows users to plot against the mouse genome
