# RegNetDriver
RegNetDriver is a computational approach to identify regulatory drivers of tumorigenesis using combined effects of coding and non-coding single nucleotide variants, structural variants (SVs) and DNA methylation changes in the DNase I hypersensitivity based regulatory network
We integrated whole-genome sequencing (WGS), RNA-Seq and DNA methylation data from primary prostate tumor samples with functional genomics data from the ENCODE and Roadmap Epigenomics projects. The three step computational model involves: (a) construction of prostate regulatory network using DNase I hypersensitivity data, (b) identification of significantly mutated, rearranged and differentially methylated coding and non-coding regulatory regions in prostate cancer and, (c) interpretation of the effects of these alterations on prostate regulatory network. We have provided scripts to construct tissue-specific regulatory network.

(A) DHS_network: folder contains source codes for tissue-specific regulatory network. It used DHS data and ENCODE promoter and enhancer annotations for creating tissue-specific DHS-based regulatory network

(B) FSig-SNV (Functionally significant single nucleotide variant) folder contains source codes to find genes with significanlty mutated coding and non-coding (promoter and enhancer) regions . It analyzes the somatic mutations in both coding and non-coding regulatory regions to identify elements that show more recurrent (present in multiple samples) and more functional mutations than expected randomly. For functional annotation, the method uses FunSeq2 to annotate and calculate functional score for each variant. Output of FSig-SNV is a list of significantly mutated coding and non-coding elements that show higher than expected frequency of functional mutations across multiple tumor samples.

(C) FSig-SV (Functionally significant structural variants) folder contains source codes to find genes with significantly rearranged coding and non-coding (promoter and enhancer) regions. FSig-SV method identifies coding and non-coding elements significantly affected by deletion, insertion, duplication, inversion and translocation events. Output of the method is a list of coding and non-coding elements that are rearranged in more samples than expected randomly. 
--To identify differentially methylated promoters and enhancer using HM450K array data, we use ELMER package (https://www.bioconductor.org/packages/release/bioc/html/ELMER.html)


How to run RegNetDriver:
COMPUTATIONAL PIPELINE FOR IDENTIFYING GENETIC AND EPIGENETIC ALTERATIONS IN THE CODING AND NON-CODING REGULATORY REGIONS OF THE TUMOR GENOME AND STUDY THEIR EFFECTS ON TISSUE_SPECIFIC REGULATORY NETWORK                                           

Steps involved:

1) Creation of Tissue-Specific Regulatory Network:

	Start by uncompressing the tissue_specific_network.tar.gz
	tar zxvf tissue_specific_network.tar.gz


	tissue_specific_network directory contains following sub-directories and files

	a) code: this folder contains following programs and scripts. spar.txt 
	PIQ_merge.pl :generates merged_files folder for PIQ output. It take PIQ output, which contains information about enrichment of TF motifs using  data.We used DHS data for prostate cell from ENCODE (https://www.encodeproject.org/experiments/ENCSR000EPU/) and encode TF motifs in JASPAR format (data/encode-motifs-Jaspar.txt ). The output of PIQcontains *-calls.csv and *.RC-calls.csv files. This program combines forward and reverse strand regions of motif enrichment into one file. 
	  
	PIQ_merged_files.sh : this script runs PIQ_merge.pl and takes directory containing PIQ output as input parameter.

	STEP2_network.sh : this scipt intersects PIQ motif predictions with tissue specific active promoters and enhancer. 

	STEP3_network.sh : this script merges multiple motifs from same family in one file. For some factors we have multiple motifs such as NR3C1_known1, NR3C1_disc1, NR3C1_known2 etc..

	STEP4_network.sh : this script generates final network. Final network is saved in NETWORK_OUTPUT directory and file TF_TARGET_EDGES.txt will contains TF_Target gene interactions, TF_enhancer_edges.txt contains TF and enhancer edges, TF_promoter_edges.txt contains TF and promoter edges.

	b) data : this folder contains input files for creating network

	c) merged_files : this folder contains merged motif files from PIQ output for prostate. This is created using PIQ_merged_files.sh script. See this script for more details.

	d) NETWORK_OUTPUT: this folder contains final output of NETWORK_DHS_SCRIPT.sh. 
	TF_TARGET_EDGES.txt file has edges between transcription factor (TF) nodes and target genes (TARGET_GENES)
	TF_enhancer_edges.txt file contains edges between TF and enhancer region
	TF_promoter_edges.txt file contains edges betweeb TF and promoter region
	TF_hubs file contains names of TF hubs and their outdegree value. Top 25% of highest outdegree nodes are selected as TF hubs.

	e) example this folder contains DHS based prostate regulatory network. Generated by running following commands:

	TO CREATE PROSTATE DHS NETWORK RUN THIS COMMAND
	sh NETWORK_DHS_SCRIPT.sh PrEC-DS12088.peaks.fdr0.01.hg19.bed merged_files

2) Find genes with significantly mutated coding and non-coding regions (FSig-SNV):

	Start by downloading FSig-SNV R package from github. See README.md for installation instructions. FSig-SNV requires annotated mutations from FunSeq2. In case of prostate tumor, we have SNVs calls from WGS data of 188 primary prostate adenocarcinoma. These 188 primary prostate sample includes 124 PRAD-CA samples from ICGC (https://dcc.icgc.org/projects/PRAD-CA), 57 samples from the work of Baca et al (PMID: 23622249)and 7 samples from Berger et al. (PMID: 21307934)

	We annotated these SNVs using FunSeq2 (https://www.ncbi.nlm.nih.gov/pubmed/25273974). We used original FunSeq2 pre-built data context available at https://github.com/khuranalab/FunSeq2_DC/tree/master/data_context

	The annotated SNVs from FunSeq2 in vcf format (Output.vcf) was used as input for FSig-SNV run. See README.md for instructions to add path variables.

	Ouptut of FSig-SNV is store in user specified location. In case of prostate tumor /path/to/Prostate/ folder will contain 'result' directory which contains 'CDS' 'promoter' 'enhancer' 'lincRNA' sub-directories. 'CDS' contains output of FSig-SNV, where Prostate_outputDf_CDS_allSamples_1000.txt file contains p-values for each genes. Here 1000 refers to number of iterations given as input by user. These pvalues can be used to create qq-plots. Similarly other folders also contains files with pvalues for genes.


3) Find genes with significantly rearranged coding and non-coding regions (FSig-SV):
	Start by downloading FSig-SV R package from github. See README.md for installation instructions. FSig-SV requires SV file in bed format. See example.bed inside FSig-Sv/input/ directory. Inputs for running FSig-SV for analyzing prostate tumor: 
Somatic Structural Variants ICGC SVs (https://dcc.icgc.org/projects/PRAD-CA), Baca et al (PMID: 23622249) and Berger et al (PMID: 21307934)

4) Find genes with significanlty differentially methylated promoter and enhancer regions (ELMER):
	We have used ELMER package. Output of ELMER is a list of genes with significantly hyper-methylated and hypo-methylated enhancers and promoters. We have enhancer definitions as used in running FSig-SV and FSig-SNV http://khuranalab.med.cornell.edu/FunSeq_data/FunSeq2_DC2/data/drm.gene.bed. 

5) Output from the above 4 steps was used to identify TF hubs with significant genetic and epigenetic alterations in coding and non-coding regulatory regions. 
	For prostate tumor, we observed 6 TFhubs (ERG, TP53, ERF, SPI1, CREB3L1 and POU2F2) significanlty altered by SVs and 3 TFhubs (TFAP2A, TFAP2C and NR3C1) significantly altered by methylation. USer can use output of steps 2,3 and 4 to find TF hubs (identified in step1) with significant alterations.


 
For any questions, comments and suggestion please email ekk2003@med.cornell.edu or prd2007@med.cornell.edu
Copyright © 2016 Ekta Khurana Lab, WCMC
