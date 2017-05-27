######################################################
### Merge multiple TF motifs into a single file.
######################################################

rm   motif_ENHANCERS/*_all   motif_PROMOTERS/*_all >/dev/null

for i in `cat data/TFgene_list`
do
	  echo -e  "Merging multiple $i motifs in $i_all file" >>logfile
	  name=$i"_"
	  cd motif_ENHANCERS/ 
		
	  cat -- "$name"*".bed" | awk '{if($10>=0.7){print }}' | sort -u -k1,1 -k2,2 -k3,3 -k4,4 -k10,10 > $i"_all" 
	  cd ..

	  
	  cd motif_PROMOTERS/
	  cat -- "$name"*".bed" | awk '{if($10>=0.7){print }}' | sort -u -k1,1 -k2,2 -k3,3 -k4,4 -k10,10 > $i"_all" 
	  cd .. 
done
echo "********* multiple motif into one TF over **********"
