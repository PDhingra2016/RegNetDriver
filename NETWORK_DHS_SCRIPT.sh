
#####################################################################################################################################

# THIS IS MAIN SCRIPT FOR GENERATING TISSUE SPECIFIC REGULATORY NETWORK      ########################################################
# MODEL EXAMPLE IS PROSTATE REGULATORY NETWORK 				     ########################################################
# INPUT1: IS DHS DATA FOR TISSUE. WE ARE USING PROSTATE DHS DATA FROM ENCODE ########################################################
# INPUT2: OUTPUT OF PIQ SOFTWARE					     ########################################################

#####################################################################################################################################

#STEP1: Find active promoter and enhancer regions in tissue and merged PIQ results files

if [ "$#" -ne 2 ];then echo -e "sh NETWORK_DHS_SCRIPT.sh <DHS_file> <merged_files>\n"
else 
	echo -e "Finding active promoters and enhancers in prostate network\n" >logfile

	DHS_file="data/"$1
	echo $DHS_file
	dir_PIQ_results=$2
	promoter_file="data/gencode.v16.promoter.bed"
	enhancer_file="data/enhancer.bed"
	home=`pwd`

# FINDING ACTIVE PROMOTERS
	intersectBed -a data/gencode_promoter_uniq -b $DHS_file -wa -wb >tmp
	awk -F ";" '$1=$1' OFS="\t" tmp |awk -F "\t" '{OFS="\t";print $1,$2,$3,$4}'|sort|uniq >"active_promoters"

# FINDING ACTIVE ENHANCERS
	awk -F "\t" '{print $1"\t"$2"\t"$3}' $enhancer_file |sort|uniq >sort_uniq_enhancer.bed
	sortBed -i sort_uniq_enhancer.bed | mergeBed >merged_sort_uniq_enhancer.bed
	intersectBed -a merged_sort_uniq_enhancer.bed -b $DHS_file -wa -wb >tmp
	awk -F "\t" '{print $1"\t"$2"\t"$3}' tmp |sort|uniq > Active_enhancers_Prostate
	intersectBed -a Active_enhancers_Prostate -b $enhancer_file -wa -wb | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$7}'| sort| uniq > "active_enhancers"
	echo -e "Active promoters and active enhancers are stored in active_promoters and active_enhancers file respectively\n" >>logfile

# MERGE PIQ OUTPUT FILES REVERSE AND FORWARD STRAND IN ONE FILE
# THIS FILE HAS INFORMATION ABOUT MOTIF ENRCIHMENT IN DHS REGIONS OF TISSUE: chr, start,end,score,purity and strand information.
# OUTPUT OF THIS PROGRAM IS STORED IN merged_files DIRECTORY
# WE HAVE ALREADY CREATED merged_files FOLDER REQUIRED FOR GENERATING PROSTATE REGULATORY NETWORK
# USER CAN RUN PIQ_merged_files.sh SCRIPT INSIDE CODE FOLDER TO CREATE merged_files FOLDER

	ls $dir_PIQ_results/* >merged_file_list 
	
	echo -e "Merged files are store in $home/merged_files folder\n" >>logfile

# STEP2: Identify motifs enriched in active promoters and active enhancer regions
	sh code/STEP2_network.sh

#STEP3: Merge multiple motifs TF motifs into a single file.  This is done to make TF-gene connections
	sh code/STEP3_network.sh

#STEP4: Generate TF-target gene edges
	sh code/STEP4_network.sh

#STEP5: Calculate degree distribution for the network
	cd NETWORK_OUTPUT
	Rscript $home/code/degree_script.R TF_TARGET_EDGES.txt 
	echo -e "Outdegree and Indegree values for TF_TARGET_EDGES.txt are stored in NETWORK_OUTPUT folder\n" >>../logfile

#STEP6
#final network stats

	num_of_TF_nodes=`awk '{print $1}' TF_TARGET_EDGES.txt |sort|uniq|wc -l`
	num_of_edges=`wc -l TF_TARGET_EDGES.txt`
	num_of_prostate_promoters=`wc -l ../active_promoters`
	num_of_prostate_enhancers=`wc -l ../active_enhancers`
	echo -e "NETWORK STATS\n" >>../logfile
	echo -e "NUMBER OF TF NODES=\t$num_of_TF_nodes" >>../logfile
	echo -e "NUMBER OF EDGES=\t$num_of_edges" >>../logfile
	echo -e "NUMBER oF PROSTATE PROMOTERS=\t$num_of_prostate_promoters" >>../logfile
	echo -e "NUMBER OF PROSTATE ENHANCERS=\t$num_of_prostate_enhancers" >>../logfile
	cd $home
	echo -e "NETWORK IS READY"

#STEP7: cleaning
	rm tmp  
	rm tmp sort_uniq_enhancer.bed merged_sort_uniq_enhancer.bed Active_enhancers_Prostate
	mkdir network_files
	mv -rf active_promoters active_enhancer motif_PROMOTERS motif_ENHANCERS merged_file_list network_files
	
fi
