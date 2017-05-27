######### this script uses PIQ output and merged reverse and forward strand  *.RC-calls.csv and *-calls.csv files into one file #####
#####################################################################################################################################
if [ "$#" -ne 1 ];then echo -e "sh PIQ_merged_files.sh  <PIQ_OUTPUT_DIR>\n"
else 
	dir_PIQ_results=$1
	home=`pwd`
	cd $dir_PIQ_results
	ls *RC*.csv >"/"$home"/temp"
	cat temp | grep ".RC" |cut -d "." -f1|sort|uniq >"/"$home"/list"
	cd $home
	mkdir merged_files			 #store all merged files in this folder
	perl PIQ_merge.pl list $dir_PIQ_results
	ls merged_files/* >merged_file_list
fi
