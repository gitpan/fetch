#!/usr/bin/perl
my @in=@ARGV;
my $where_index, $question, $i, $s,$t,$fasta,$index_file,
   $all, $target, $matched, $seq, $gdf, $list, $count, $create;
if(@ARGV < 1){
  &HELP; sub HELP{
	   print "\n Usage: $0 [-cafghl] <any swissprot name entry>\n";
	   print "   -f is for FASTA output format only\n";
	   print "   -a is for ALL matched seq. \(same as using glob=\> *YEAST\)\n";
	   print "   -c is for Creating seq.idx file\n";
	   print "   -h is for HELP!\n";
	   print "   -g is for GDF file format output\n";
	   print "   -l is for list of match entries\(in 1 column\)\n";
	   print "\n\n  $0 utilizes 2 ENV setting if exist: SWDIR or SWISS and SWINDEX";
	   print "\n   You have to set ENV vars, SWDIR or SWISS to seq.dat path\n";
	   print "\n   In  csh: \'setenv SWDIR /your/DB/swiss\'";
	   print "\n   In bash: \'export SWDIR=/your/DB/swiss\'\n";
	   print "\n   Option: you can also set ENV var, SWINDEX to seq.idx path\n";
	   print "   You can put absolute path for seq.dat and(or) seq.idx at prompt\n\n\n";
	   print chr(7);} }else{ for($t=0; $t<@in; $t++){if($in[$t]=~/^\-h$/){ &HELP; exit;}if($in[$t]=~/^\-c$/i){$create=1; splice(@in, $t, 1); $t--;
		   print "\n You should have provided \/path\/seq.dat file with this option\n"; next; }
		if($in[$t]=~/^\-af$/){ $fasta=$all=1; splice(@in, $t, 1); $t--; next;}
		if($in[$t]=~/^\-ag$/){ $gdf=$all=1; splice(@in, $t, 1); $t--; next;}
		if($in[$t]=~/^\-g$/){    $gdf=1; splice(@in, $t, 1); $t--; next}
		if($in[$t]=~/^\-f$/i){   $fasta=1; splice(@in, $t, 1); $t--; next}
		if($in[$t]=~/^\-a$/i){   $all=1;   splice(@in, $t, 1); $t--; next}
		if($in[$t]=~/^\-l$/i){   $list=$all=1;   splice(@in, $t, 1); $t--; next}
		if( ($in[$t]=~/seq\.dat/)&&(-f $in[$t])){ $DB=$in[$t];  splice(@in, $t, 1); $t-- }
		if( ($in[$t]=~/seq\.idx/)&&(-e $in[$t])){ $index_file=$in[$t];	splice(@in, $t, 1); $t--;}}
	 if(($create==1)&&(defined($DB)) ){ goto CREATE; } elsif($create==1){ print "\n You must give seq.dat file with path as an arg\n"; exit }   }
   unless(defined($DB)){if(-e "$ENV{'SWDIR'}seq.dat" ){ $DB="$ENV{'SWDIR'}seq.dat";
	  }elsif( -f "$ENV{'SWISS'}seq.dat" ){$DB="$ENV{'SWISS'}seq.dat";
	  }elsif( -e 'seq.dat'){$DB="seq.dat"; }elsif( -f "$ENV{'swiss'}seq.dat"){$DB="$ENV{'swiss'}seq.dat";
	  }elsif(-f "ENV{'HOME'}seq.dat"){$DB="ENV{'HOME'}seq.dat"; }elsif(-f "ENV{'SWDIR'}\/seq.dat"){$DB="ENV{'SWDIR'}\/seq.dat";
	  }else{ ASK: print "\n Where is your swissprot seq.dat file?\n";
		 $swiss=<STDIN>; chomp($swiss);if(-e "$swiss"){goto SW_INDEX;
		 }else{	goto ASK}}} unless(defined($index_file)){ SW_INDEX:
	  if(-e "$ENV{'SWDIR'}seq.idx" ){ $index_file="$ENV{'SWDIR'}seq.idx";
	  }elsif( -f "$ENV{'SWISS'}seq.idx" ){ $index_file="$ENV{'SWISS'}seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}" ){$index_file="$ENV{'SWISS'}";
	  }elsif( -e 'seq.idx'){ $index_file="seq.idx";
	  }elsif( -f "$ENV{'swiss'}seq.idx"){$index_file= "$ENV{'swiss'}seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}seq.idx"){$index_file= "$ENV{'SWINDEX'}seq.idx";
	  }elsif( -f "$ENV{'HOME'}seq.idx"){$index_file= "$ENV{'swiss'}seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}seq.idx"){$index_file="$ENV{'SWINDEX'}\/seq.idx";
	  }elsif( -f "$ENV{'swindex'}seq.idx"){$index_file="$ENV{'swindex'}seq.idx";
	  }elsif(defined($DB)){ print "\n Your swissprot is in $DB, but no seq.idx file from it.\n";
		 W: print "\n Where is your seq.idx \(type path and filename\)?\n";
		 $question++; $where_index=<STDIN>; chomp($where_index); if(-f $where_index){
			open(TMP, "$where_index"); while(<TMP>){ if($_=~/^ *\w+ +\d+/){ $index_file=$where_index;
				   print "\n Your index file seems to be right \($index_file\) \n";
				   goto MAIN_SEARCH; }elsif($count > 4){ # read at least 4 lines and see if they are index
				   print "\n $where_index doesn't seem to be index file\n";
				   print "\n Terminate(t) or go on (g) trying\n";
				   $try=getc; if($try=~/t/i){  exit; } else{ goto W; }}else{ $count++;}} close TMP; }else{
			if($question > 3){ print "\n I can create the index in pwd and run $0 and \n";
			   print "\n you can copy seq.idx into your swissprot dir later\n";
			   goto CREATE;	}goto W; } CREATE: print "\n Can I create seq.idx in pwd? (y+return or return)\n";
		 $yes_no=getc;if($yes_no=~/y/i){print "\n seq.dat is being created....\(1 minute in my Linux\)\n";
			open(DB, "$DB");open(IDX, ">seq.idx");	while(<DB>){
			  if(/^ID[\t ]+(\w+) +/){ $index{$1}=tell(DB);
				 print IDX "\n$1 $index{$1}";}}	close(DB, IDX);
			if(-s "seq.idx"){ print "\nGood. seq\.idx file is created.";
			    print "\n Copy seq.idx at SWISSPROT dir or you can set\n";
				print "absolute path ENV var \'SWINDEX\' to your seq.idx path\n";
				print "e.g. #bash\> export SWINDEX=\/DB\/Swiss\/seq.idx\n\n";
				if($create==1){ exit}}else{ print "\n Creation of seq.dat seems to have gone wrong"; }
		 }else{	exit }}}  MAIN_SEARCH:  for($i=0; $i<@in; $i++){
	  my @possible, @pos, %possible;  my $target=$in[$i];
	  if($target=~/\*/){$target=~s/\*/\\\w\{0,4\}/; # to handle glob input
		 $all=1; } open(INDEX, "$index_file");
	  while(<INDEX>){if( /(\w*$target\w*) +(\d+)/ ){
		   $possible{$1}=$2;} } close INDEX;
	  @poss = keys %possible; if( (@poss >1)&&($all !=1)){
		 print "\n @poss","\n";	print chr(7); print "\n There are more than a few seqs for $in[$i]";
		 print "\n be more specific! OR use -a option for all matched\n\n";
	  }elsif($all !=1){ print "\n"; open (DB, "$DB"); A:for($p=0; $p < @poss; $p++){
		   if($poss[$p]=~/\w*$target\w*/){ $matched=$possible{$poss[$p]};
			 seek(DB, ($matched-51), 0); while(<DB>){ if($gdf==1){
			      if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){ printf ("%-24s %-3d %-7d %-14s %4s\n", "$poss[$p]\/1\-$1", 1, $1, $poss[$p], '0.0');
					 next A;  } }elsif((/^\/\/$/)&&($fasta==1)){ $seq=~s/ //g; print "\>$poss[$p]\n$seq\n";
				  $seq=''; next A; }elsif((/^\/\/$/) && ($fasta !=1)){  print "\/\/\n"; next A; }elsif(($fasta==1)&&(/^[\t ]+\w+/)){ $seq.=$_; next ; }elsif($list==1){ if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         print "$poss[$p]\n"; next A; } }elsif($fasta !=1){ print }}}} close(DB);
	  }elsif($all==1){ print "\n"; open (DB, "$DB"); A:for($p=0; $p < @poss; $p++){
		   if($poss[$p]=~/\w*$target\w*/){$matched=$possible{$poss[$p]};
			 seek(DB, ($matched-51), 0); while(<DB>){if($gdf==1){if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         printf ("%-24s %-3d %-7d %-14s %4s\n", "$poss[$p]\/1\-$1", 1, $1, $poss[$p], '0.0');
					 next A;} }elsif((/^\/\/$/)&&($fasta==1)){$seq=~s/ //g;
				  print "\>$poss[$p]\n$seq\n"; $seq=''; next A; }elsif((/^\/\/$/) && ($fasta !=1)){
				  print "\/\/\n"; next A; }elsif(($fasta==1)&&(/^[\t ]+\w+/)){
				  $seq.=$_;  next; }elsif($list==1){ if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         printf "$poss[$p]\n"; next A; } }elsif($fasta !=1){ print ;}}}} close(DB);} }

