# STEP1: reads list of files created by PIQ and use *-calls.csv, *.RC-calls.csv and *-calls.all.bed , *.RC-calls.all.bed files generated from PIQ run to generate one single motif file with merged output motif enrichment results for forward  and reverse strands.
#the output from this step will be used to intersect active promoter and enhancer regions.
# Input: list
# Output: <motifname>-calls.merged.filtered.bed
# Execute: perl PIQ_merge.pl list
#
#created by : Priyanka Dhingra Aug 23, 2016
#
use Data::Dumper;
use File::Basename;
if (@ARGV!=2) {
    print "STEP2_merge.pl [File containing lists of base names of PIQ] [PIQ directory]";
}
else{
    my $FP;
    my %motif;
    my @keys;my $x=0;
    my @keys=(		  "TFname",
                          "PIQ_TF",
                          "encode",
			  "motif",
			  "suffix",
             );
    open($FP,"motif_names_file.txt") or die("Cannot open motif_names_files.txt\n");
    while(<$FP>){
	my @arr_name;
	chomp($_);
	@arr_name=split(' ',$_);
	$motif{TFname}[$x]=$arr_name[0];	
	$motif{PIQ_TF}[$x]=$arr_name[1];	
	$motif{encode}[$x]=$arr_name[2];	
	$motif{motif}[$x]=$arr_name[3];	
	$motif{suffix}[$x]=$arr_name[4];
	$x++;	
    }
    close($FP);

    my $root_list_fh;
    open($root_list_fh, '<', $ARGV[0]) or die ("Cannot open $ARGV[0]!");
    while (my $root_line=<$root_list_fh>) {

        chomp($root_line);
	my($filename, $dirs, $suffix) = fileparse($root_line);

	my $dir=$ARGV[1];
	my @line_array=$filename=~m/([A-Za-z0-9]+)/g;
        my $motif_name=$line_array[1];

	for($i=0;$i<scalar(@{$motif{TFname}});$i++){
		if($motif{PIQ_TF}[$i] eq $motif_name) {$dirs=$motif{TFname}[$i]; last;}
		}

        my $forward_csv=$dir."/".$root_line."-calls.csv";
        my $forward_bed=$dir."/".$root_line."-calls.all.bed";
        my $reverse_csv=$dir."/".$root_line.".RC-calls.csv";
        my $reverse_bed=$dir."/".$root_line.".RC-calls.all.bed";
        my $out_bed="merged_files/".$dirs."_".$motif_name."-calls.merged.filtered.bed";
	my $fh; 
	print $motif_name." ".$out_bed."\n";
	
	%forward_info=read_file($forward_csv);
	%reverse_info=read_file($reverse_csv);

	my $outfile_stream;
	open($outfile_stream, '>', $out_bed) or die ("Cannot open $out_bed for writing!");
#	print $outfile_stream "track type=\"BAM\" db=\"hg19\" name=\"$motif_name\" description=\"Hits for motif $motif_name\" usescore=1\n";
	read_bed_entries($forward_bed, $outfile_stream, \%forward_info, $motif_name);
	read_bed_entries($reverse_bed, $outfile_stream, \%reverse_info, $motif_name);
	close $outfile_stream;
}
}

sub read_file{
	my $input_file=$_[0];
	my %read_file;
	my $x=0; my @keys;
	 my @keys=("chr",
                          "coord",
                          "pwm",
			  "shape",
			  "score",
			  "purity"
                         );

	my $y=0;
	open(FA,$input_file) or die("Cannot open file $input_file");
	while(my $line=<FA>)
	{
			chomp($line);
			my @contents=split(',',$line);
			tr/"//d foreach(@contents);

			$read_file{chr}[$y]=$contents[1];
			$read_file{coord}[$y]=$contents[2];
			$read_file{pwm}[$y]=$contents[3];
			$read_file{shape}[$y]=$contents[4];
			$read_file{score}[$y]=$contents[5];
			$read_file{purity}[$y]=$contents[6];
			$y++;
	}
	close(FH);
	#print Dumper(\%read_file);
	return(%read_file);

}
sub read_bed_entries{
    my $infile_name=$_[0];
    my $out_bed_object=$_[1];
    my %info_hash=%{$_[2]};
    my $motif_id=$_[3];
    my $infile_fh;
    open($infile_fh, '<', $infile_name) or die ("Cannot opne $infile_name for reading!");
    my $first_line=<$infile_fh>;
    my $counter=0;
    while( my $one_feature =<$infile_fh>){
	chomp($one_feature);
	my @feature_array=split '\t', $one_feature;
	for(my $counter=0;$counter<scalar(@{$info_hash{chr}});$counter++)
	{	   	
           if (($feature_array[1] eq $info_hash{coord}[$counter])&& ($feature_array[0] eq $info_hash{chr}[$counter])) 
 	  {
	    printf $out_bed_object "$feature_array[0]\t$feature_array[1]\t$feature_array[2]\t$motif_id\t$info_hash{score}[$counter]\t$info_hash{purity}[$counter]\t$feature_array[5]\n";
	    last;
           }
	}
            
        
    }
    close($inflie_fh);
}
