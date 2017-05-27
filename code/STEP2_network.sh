##################################################################################################
##### INTERSECT PIQ MOTIF PREDICTIONS IN DHS REGIONS WITH ACTIVE PROMOTERS & ENHANCERS  ##########
##################################################################################################
 
rm -rf motif_PROMOTERS motif_ENHANCERS >/dev/null

mkdir motif_PROMOTERS
mkdir motif_ENHANCERS
echo -e  " Identifying motifs enriched in active_promoters and active_enhancers" >>logfile
for i in `cat merged_file_list`
do

	 name=`basename $i -calls.merged.filtered.bed`	
	 intersectBed  -a active_promoters -b $i -wa -wb  >"motif_PROMOTERS/"$name"-intersect.bed"
	 intersectBed  -a active_enhancers -b $i -wa -wb  >"motif_ENHANCERS/"$name"-intersect.bed"
done	
echo -e  "Motifs enriched in active promoters and active enhancers are stored in motif_PROMOTERS and motif_ENHANCERS folder\n" >>logfile
