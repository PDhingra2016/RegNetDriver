################################################################
#####   SCRIPT TO MAKE FINAL NETWORK (TF-target gene edges)
###############################################################

rm -rf NETWORK_OUTPUT
mkdir NETWORK_OUTPUT

dir="NETWORK_OUTPUT"

rm ENH_DHS_NETWORK  PRO_DHS_NETWORK

for i in `cat data/TFgene_list`
do 
	
	awk  -v i=$i -F "\t" '{OFS=FS;  print $1,$2,$3,$4,i,$10,"Enhancer"}' motif_ENHANCERS/$i"_all" >> $dir"/ENH_DHS_NETWORK"	
	awk -v i=$i  -F "\t" '{OFS=FS;  print $1,$2,$3,$4,i,$10,"Promoter"}' motif_PROMOTERS/$i"_all" >> $dir"/PRO_DHS_NETWORK"	
done


cd $dir

sort -u  ENH_DHS_NETWORK >unique_ENH_DHS_NETWORK
sort -u  PRO_DHS_NETWORK >unique_PRO_DHS_NETWORK
echo "TF\tchr_pro\tstart_pro\tend_pro" >TF_promoter_edges.txt
echo "TF\tchr_enh\tstart_enh\tend_enh" >TF_enhancer_edges.txt
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5}' unique_ENH_DHS_NETWORK |sort| uniq >>TF_enhancer_edges.txt
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$5}' unique_PRO_DHS_NETWORK |sort| uniq >>TF_promoter_edges.txt

cat unique_ENH_DHS_NETWORK unique_PRO_DHS_NETWORK >FINAL_ENH_PRO_DHS_NETWORK
echo -e "TF - target gene edges with promoter and enhancer information are stored in FINAL_ENH_PRO_DHS_NETWORK file\n" >>../logfile
echo  "TF\tTARGET_GENE" >TF_TARGET_EDGES.txt
awk '{print $5"\t"$4}' FINAL_ENH_PRO_DHS_NETWORK | sort |uniq >>TF_TARGET_EDGES.txt
echo -e "FINAL PROSTATE NETWORK WITH TF-TARGET GENE EDGES IS STORED IN TF_TARGET_EDGES.txt\n" >>../logfile

cd ../
