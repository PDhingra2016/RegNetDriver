#!/bin/bash
#ICGC_data=./ICGC_structural_somatic_mutation.PRAD-CA.tsv
if [ "$ICGC_data" == "" ]; then
        echo "ICGC data not specified."
else
awk '{ OFS="\t"
#        if($9 ~/del/){ print "chr"$10,$11,$14,"del",$9,$5}
#        else if ($9 ~/inv/){ print "chr"$10,$11,$14,"inv",$9,$5}
#        else if ($9 ~/dup/){ print "chr"$10,$11,$14,"dup",$9,$5}
#        else if (($9 ~/tra/) && ($10 == $13)){ print "chr"$10,$11,$14,"tra1",$9,$5}
#        else if (($9 ~/tra/) && ($10 != $13)){ print "chr"$10,$11-10000,$11+10000,"tra2",$9,$5"\nchr"$13,$14-10000,$14+10000,"tra2",$9,$5}
#}' ICGC_structural_somatic_mutation.PRAD-CA.tsv > ./input/ICGC_SV_data.bed
fi
