#!/usr/bin/perl
# Last Update by ./update_subroutines.pl: Thu Jan  2 03:25:32 GMT 1997
#______________________________________________________________
# Title     : fetch_seq.pl
# Usage     : fetch.pl SWISS_ENTRY
# Function  : fetches swissprot entry or fasta format seq with
#             given seq name(like  SAA_HORSE, SA*HORSE, SAA,..)
#             you can give multi files(SAA*, SAU*) at the same
#             time. This uses ENV setting of 'SWDIR'.
#             If you already know where you put seq.dat you can
#             specify in the prompt.  jong@mrc-lmb.cam.ac.uk
# Example   :  fetch_seq.pl YAPT
# Keywords  : fetch_swissprot_sequence, fetch_sequence,
#             find_swiss_sequence, find_sequence
# Options   : _  for debugging.
#             #  for debugging.
#             -f for fasta format file output
#             -a is for ALL matched seq. (same as using glob=> *YEAST)
#             -c is for Creating seq.idx file
#             -h is for HELP!
#             -g is for GDF file format output
#             -l is for list of match entries(in 1 column)
#             -s is for species option (input name mst be species (YEAST, RAT, HUMAN..)
#             n= is for Number of seq you want to get from swissprot
#             s= is for Size limit. Min seq size in swiss, s=10  -> minimum 11 aa seq.
#             S= is for Size limit. Max seq size in swiss, s=1000 -> get less than 1000
#
# Returns   : STDOUT
# Argument  : swissprot seqname
# Version   : 1.6
#--------------------------------------------------------------


&fetch_seq(@ARGV);


#______________________________________________________________
# Title     : fetch_seq
# Usage     : &fetch_seq(@ARGV);
# Function  : fetches swissprot entry or fasta format seq with
#             given seq name(like  SAA_HORSE, SA*HORSE, SAA,..)
#             you can give multi files(SAA*, SAU*) at the same
#             time. This uses ENV setting of 'SWDIR'
# Example   : &fetch_swiss_seq(@ARGV);
# Keywords  : fetch_swissprot_sequence, fetch_sequence,
#             find_swiss_sequence, find_sequence
# Options   : _  for debugging.
#             #  for debugging.
#             -f for fasta format file output
#             -a is for ALL matched seq. (same as using glob=> *YEAST)
#             -c is for Creating seq.idx file
#             -h is for HELP!
#             -g is for GDF file format output
#             -l is for list of match entries(in 1 column)
#             -s is for species option (input name mst be species (YEAST, RAT, HUMAN..)
#             n= is for Number of seq you want to get from swissprot
#             s= is for Size limit. Min seq size in swiss, s=10  -> minimum 11 aa seq.
#             S= is for Size limit. Max seq size in swiss, s=1000 -> get less than 1000
#
# Argument  : swissprot seqname
# Version   : 1.6
#--------------------------------------------------------------
sub fetch_seq{
   my @in=@_;
   my $FASTA_index, $FASTA, $where_index, %index, $question, $i,
	  $s,$t,$fasta,$index_file, $all, $species,$target, $matched, $seq, $gdf, $list, $count, $create;
   my $SEQ_size_max=100000000;

   if(@_ < 1){	  &HELP_fetch_seq;
   }else{
	 F: for($t=0; $t<@in; $t++){ #'''''''''''' PROMPT ARGV processing ''''''''''''''''''
		if($in[$t]=~/^\-c$/i){
		   $create=1; splice(@in, $t, 1); $t--;
		   print "\n You should provide database\(e.g, seq.dat\) file with this opt, I guess you did\n";
		   print "\n If you wanted to make an index with any fasta db, you also have to\n";
		   print "  give the file name. e.g:\n     $0 -c /DB/swiss/seq.dat\n";
		   print "  or $0 -c my_db.fa\n\n";
		   next; }
		if($in[$t]=~/^\-af$/){ $fasta=$all=1; splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-afs$/){ $species=$fasta=$all=1; splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-ag$/){ $gdf=$all=1; splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-g$/){    $gdf=1; splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-f$/i){   $fasta=1; splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-a$/i){   $all=1;   splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-l$/i){   $list=$all=1;   splice(@in, $t, 1); $t--; next; }
		if($in[$t]=~/^\-s$/i){  $species=$all=1; splice(@in, $t, 1); $t--; next; }
		if( ($in[$t]=~/seq\.dat/)&&(-f $in[$t])){ ## if the path for swiss prot is given
		   $DB=$in[$t];  splice(@in, $t, 1); $t--; next;        }
		if( ($in[$t]=~/seq\.idx/)&&(-e $in[$t])){ ## if the path for swiss index is given
		   $index_file=$in[$t];	splice(@in, $t, 1); $t--; next;	}

		#''''''' SWiss prompt input file check ''''''''''''''''''
		if( -f $in[$t]){
		   open(TEMP, "$in[$t]");
		   while(<TEMP>){
			 if(/^ID[\t ]+\w+/){$DB=$in[$t]; splice(@in, $t, 1);$t--;next F;}}
		   close TEMP;
		}

		#'''''''' FASTA prmpt input file check '''''''''''''''''''
		if( (-f $in[$t]) && !(defined($FASTA))){
		   open(TEMP, "$in[$t]");
		   while(<TEMP>){
			 if(/^\> {0,4}\S+/){$FASTA=$in[$t]; $fasta=1;
			 if(-s "$FASTA\.idx"){ $FASTA_index="$FASTA\.idx"; }
		     splice(@in, $t, 1);$t--;next F;}}
		   close TEMP;
		}

		#'''''''' INDEX file automatic check ''''''''''''''''''
		if( -f $in[$t]){
		   open(TEMP2, "$in[$t]");
		   my $first_pos, $Count, @splited;
		   while(<TEMP2>){
			 $Count++;
			 if( $Count>3 ){
				if(/^ {0,2}\S+ +(\d+)/){
				   if(defined($first_pos) && ($1-$first_pos ) > 1000 ){
					  $index_file=$in[$t];
					  splice(@in, $t, 1);$t--;next F;
				   }elsif( defined($first_pos) && ($1-$first_pos)<1000 ){
					  $FASTA_index=$in[$t]; $fasta=1;
					  if($FASTA_index=~/^(\S+)\.\w+$/){
					     if(-s $1){ $FASTA= $1; }
					  }
					  splice(@in, $t, 1);$t--;next F;
				   }
				   $first_pos=$1;
				}
			 }
		   }
		   close TEMP2;
		}
		if($in[$t]=~/^\-h$/i){ &HELP_fetch_seq; exit;}
		if($in[$t]=~/^n=(\d+)$/i){ $SEQ_num_to_fetch=$1;
		   splice(@in, $t, 1);$t--;next F;}
		if($in[$t]=~/^s=(\d+)$/){ $SEQ_size_min=$1; $fasta=1;
		   splice(@in, $t, 1);$t--;next F;}
		if($in[$t]=~/^S=(\d+)$/){ $SEQ_size_max=$1; $fasta=1;
		   splice(@in, $t, 1);$t--;next F;}
	 }

	 if(($create==1)&&(defined($DB)) ){ goto CREATE; }
	 elsif(($create==1) && (defined($FASTA)) ){ goto CREATE; }
	 elsif($create==1){
	    print "\n You must give db filename (e.g. seq.dat) with path to make an index";
	    print "\n  I can handle fasta db file to make an index\n";
	    exit;
	 }
   }

   if($SEQ_size_max < $SEQ_size_min){ print "\n Seq size Max is smaller than min\n"; exit; }

   ##""""""""""""""""""""""" DB file if not defined """"""""""""""""""""""""""""""""""""""""""""
   if (!defined($DB)){
	  if((!defined($FASTA))&&($fasta==1)&&(-T "$ENV{'FASTADB'}")){
		 $FASTA=$ENV{'FASTADB'};
	  }elsif(defined($FASTA) && ($fasta==1) &&($create !=1) ){
		 goto SW_INDEX;
	  }elsif(!defined($FASTA) && (defined($FASTA_index))&& !(-e "$ENV{'FASTADB'}") ){
		 print "\n NO fasta db is defined\n";
		 goto ASK;
	  }elsif(-e "$ENV{'SWDIR'}seq.dat" ){
		 $DB="$ENV{'SWDIR'}seq.dat";
	  }elsif(-e "$ENV{'FETCHSWISS'}seq.dat" ){
		 $DB="$ENV{'FETCHSWISS'}seq.dat";
	  }elsif(-e "$ENV{'FETCHSWISS'}" ){
		 $DB="$ENV{'FETCHSWISS'}";
	  }elsif(-e "$ENV{'SWDIR'}\/seq.dat" ){
		 $DB="$ENV{'SWDIR'}\/seq.dat";
	  }elsif( -f "$ENV{'SWISS'}seq.dat" ){
		 $DB="$ENV{'SWISS'}seq.dat";
	  }elsif( -f "$ENV{'SWISS'}\/seq.dat" ){
		 $DB="$ENV{'SWISS'}\/seq.dat";
	  }elsif( -e 'seq.dat'){
		 $DB="seq.dat";
	  }elsif( -f "$ENV{'swiss'}seq.dat"){
		 $DB="$ENV{'swiss'}seq.dat";
	  }elsif(-f "ENV{'HOME'}seq.dat"){
		 $DB="ENV{'HOME'}seq.dat";
	  }elsif(-f "ENV{'SWDIR'}\/seq.dat"){
		 $DB="ENV{'SWDIR'}\/seq.dat";
	  }else{
		ASK: print "\n Where is your swissprot seq.dat(or fasta db) file?\n";
			 print "  I recommand you to set the path for them in ENV vars\n";
			 print "  e.g. export SWDIR=/DB/Swiss/  to where you put seq.dat\n";
			 print "  e.g. export FASTADB=/DB/Swiss/my_swiss.fa  for fasta database\n";
		 $swiss=<STDIN>;
		 chomp($swiss);
		 if( -f $swiss){
			open(TEMP, "$swiss");
			while(<TEMP>){
			   if(/^ID[\t ]+\w+/){ $DB=$swiss; goto SW_INDEX; }
			   elsif(/^\> {0,3}\S+/){ $FASTA=$swiss; goto SW_INDEX;}
			}
			close TEMP;
		 }else{
			goto ASK;
		 }
	  }
   }
   ##""""""""""""""""""""""""""""" INDEX file """"""""""""""""""""""""""""""""""""""""
   if( !defined($index_file)){
	  SW_INDEX:
	  if((!defined($FASTA_index))&&($fasta==1)&&(-T "$ENV{'FASTAINDEX'}")){
		 $FASTA_index=$ENV{'FASTAINDEX'};
	  }elsif(!defined($FASTA_index)&&(-T $FASTA)){
		 goto W;
	  }elsif(defined($FASTA_index)&&(-T $FASTA)){
		 goto MAIN_SEARCH;
	  }elsif(-e "$ENV{'FETCHSWISSINDEX'}seq.idx" ){
		 $index_file="$ENV{'FETCHSWISSINDEX'}seq.idx";
	  }elsif(-e "$ENV{'FETCHSWISSINDEX'}\/seq.idx" ){
		 $index_file="$ENV{'FETCHSWISSINDEX'}\/seq.idx";
	  }elsif(-e "$ENV{'SWDIR'}seq.idx" ){
		 $index_file="$ENV{'SWDIR'}seq.idx";
	  }elsif( -f "$ENV{'SWISS'}seq.idx" ){
		 $index_file="$ENV{'SWISS'}seq.idx";
	  }elsif( -f "$ENV{'SWISS'}\/seq.idx" ){
		 $index_file="$ENV{'SWISS'}\/seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}" ){
		 $index_file="$ENV{'SWINDEX'}";
	  }elsif( -e 'seq.idx'){
		 $index_file="seq.idx";
	  }elsif( -f "$ENV{'swiss'}seq.idx"){
		 $index_file= "$ENV{'swiss'}seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}seq.idx"){
		 $index_file= "$ENV{'SWINDEX'}seq.idx";
	  }elsif( -f "$ENV{'HOME'}seq.idx"){
		 $index_file= "$ENV{'HOME'}seq.idx";
	  }elsif( -f "$ENV{'SWINDEX'}seq.idx"){
		 $index_file="$ENV{'SWINDEX'}\/seq.idx";
	  }elsif( -f "$ENV{'swindex'}seq.idx"){
		 $index_file="$ENV{'swindex'}seq.idx";
	  }elsif(defined($DB)|| defined($FASTA) ){
		 print "\n Your swissprot is in $DB, but no seq.idx file for it.\n";
		 W: print "\n Where is seq.idx(or fasta idx file eg. $FASTA\.idx), type path and filename?\n";
		    print "  I recommand you to set the path for them in ENV vars later\n";
			print "  e.g. export SWINDEX=/DB/Swiss/  to where you put seq.dat index\n";
			print "  e.g. export FASTAINDEX=/DB/Swiss/my.fa.idx  for fasta db index\n";
			print "  Asking where 3 times. After, will ask creation of seq.idx or $FASTA.idx\n";
		 $question++;
		 $where_index=<STDIN>;
		 chomp($where_index);
		 if(-f $where_index){
			open(TMP, "$where_index");
		    while(<TMP>){
				if($_=~/^ {0,2}\S+ +\d+/){
				   $index_file=$where_index;
				   print "\n Your index file seems to be right \($index_file\) \n";
				   goto MAIN_SEARCH;
				}elsif($count > 4){ # read at least 4 lines and see if they are index
				   print "\n $where_index doesn't seem to be index file\n";
				   print "\n Terminate(t) or go on (g) trying\n";
				   $try=getc;
				   if($try=~/t/i){  exit; }
				   else{ goto W; }
				}else{
				   $count++;
				}
			}
			close TMP;
		 }else{
			if($question > 2){
			   print "\n I can create the index in pwd for you run $0 and \n";
			   print "\n you can copy seq.idx(or $FASTA\.idx) into your swissprot dir later\n";
			   goto CREATE;
			}
			goto W;
		 }

		 #""""""""""""""" CREATION of INDEX file """""""""""""""""""""""""""""""""""""""""""""
		 CREATE:
		 if(defined($DB)){ print "\n Can I create seq.idx in pwd? (y+return or return)\n" }
		 if(defined($FASTA)){ print "\n Can I create $FASTA\.idx in pwd? (y+return or return)\n" }
		 $yes_no=getc;
		 if($yes_no=~/y/i){
			if(defined($DB)){
			   print "\n seq.idx being created...\(1 min in my Linux\)\n";
			   open(DB, "$DB");
			   open(IDX, ">seq.idx");
			   print IDX "# swiss_index\n";
			   while(<DB>){
				 if(/^ID[\t ]+(\w+) +/){
					$index{$1}=tell(DB);
					print IDX "\n$1 $index{$1}";
				 }
			   }
			   close(DB, IDX);
			   if(-s "seq.idx"){
				   print "\nGood. seq\.idx is created.";
				   print "\n Copy seq.idx to SWISSPROT dir or you can set\n";
				   print "absolute path ENV var \'SWINDEX\' to your seq.idx path\n";
				   print "e.g. #bash\> export SWINDEX=\/DB\/Swiss\/seq.idx\n\n";
				   if($create==1){ exit;  }
			   }else{
				   print "\n Creation of seq.dat seems to have gone wrong";
			   }

			}elsif(defined($FASTA)){
			   $F_idx="$FASTA\.idx";
			   print "\n $F_idx being created...\n";
			   open(FASTADB, "$FASTA");
			   open(FASTAIDX, ">$F_idx");
			   print FASTAIDX "# fasta_index\n";
			   while(<FASTADB>){
				 if(/^\> {0,4}(\S+) */){
					$index{$1}=tell(FASTADB);
					print FASTAIDX "\n$1 $index{$1}";
				 }
			   }
			   close(FASTADB, FASTAIDX);
			   if(-s $F_idx){
				   print "\nGood! Copy $F_idx to your DB dir and set two ENV vars\n";
				   print "absolute path ENV var \'FASTADB\' to your fastadb path\n";
				   print "absolute path ENV var \'FASTAINDEX\' to your $F_idx path\n";
				   print "e.g. #bash\> export FASTADB   =\/DB\/mySwiss\/$FASTA\n";
				   print "e.g. #bash\> export FASTAINDEX=\/DB\/mySwiss\/$F_idx\n";
				   print "e.g. #tcsh\> setenv FASTADB    \/DB\/mySwiss\/$FASTA\n";
				   print "e.g. #tcsh\> setenv FASTAINDEX \/DB\/mySwiss\/$F_idx\n";
				   print "Unless, you can specify the database each time at prompt\n\n";
				   if($create==1){ exit;  }
			   }else{
				   print "\n Creation of seq.dat or $F_idx seems to have gone wrong";
			   }
			}
		 }else{
			exit;
		 }
	  }
   }

   #""""""""""""""""""""""""""" MAIN SERACH """""""""""""""""""""""""""""""""""""""""""""""
   MAIN_SEARCH:
   for($i=0; $i<@in; $i++){
	  my (@possible, @pos, %possible); my $target=$in[$i];
	  if($target=~/\*/){
		 $target=~s/\*/\\\w\{0,6\}/; # to handle glob input
		 $all=1;
	  }
	  if(defined($index_file)){
		 open(INDEX, "$index_file");
		 if($species==1){
		    while(<INDEX>){
		      if( /(\w*\_$target) +(\d+)/ ){ $possible{$1}=$2; }
		    }
		 }else{
		    while(<INDEX>){
		      if( /(\w*$target\w*) +(\d+)/ ){ $possible{$1}=$2; }
		    }
		 }
		 close INDEX;
		 goto SWISS;
	  }elsif(($fasta==1) && (defined($FASTA_index)) ){
		 open(INDEX, "$FASTA_index");
		 if($species==1){
		    while(<INDEX>){
		      if( /(\w*\_$target) +(\d+)/ ){ $possible{$1}=$2; }
		    }
		 }else{
		    while(<INDEX>){
		      if( /(\w*$target\w*) +(\d+)/ ){ $possible{$1}=$2; }
		    }
		 }
		 close INDEX;
		 goto FASTA;
	  }

	  SWISS:
	  @poss = sort keys %possible;

	  if( (@poss >1)&&($all !=1)){
		 print "\n @poss","\n";
		 print chr(7);
		 print "\n There are more than a few seqs for $in[$i]";
		 print "\n be more specific! OR use -a option for all matched\n\n";
		 exit;
	  }elsif($all !=1){
		 print "\n";
		 open (DB, "$DB");
		 if(defined($SEQ_num_to_fetch)){
			print "\n# You defined the number of sequence to fetch: $SEQ_num_to_fetch\n";
			$num_sequence=$SEQ_num_to_fetch;
		 }else{ $num_sequence=@poss; }

		 A:for($p=0; $p < $num_sequence; $p++){
		   if($poss[$p]=~/\w*$target\w*/){
			 $matched=$possible{$poss[$p]};   # %possible has the name and index num
			 seek(DB, ($matched-52), 0);
			 while(<DB>){
			   if($gdf==1){
			      if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         printf ("%-24s %-3d %-7d %-14s %4s\n", "$poss[$p]\/1\-$1", 1, $1, $poss[$p], '0.0');
					 next A;
			      }
			   }
			   elsif(/^ {0,2}\/\// and  $fasta !=1){  # !!! DO NOT put $ in /^ {0,2}\/\// as there is something
				  print "\/\/\n";
				  next A;
			   }elsif(/^ {0,2}\/\//  and  $fasta==1){ # !!! DO NOT put $ in /^ {0,2}\/\// as there is something
				  $seq=~s/ //g;
				  if( ($SEQ_size_min < length($seq))&&(length($seq) < $SEQ_size_max) ){
					 print "\>$poss[$p]\n$seq\n"; $seq=''; next A;
				  }else{  $seq=''; $num_sequence++;  next A; }
			   }elsif( $fasta==1 and /^[\t ]+\w+/){
				  $seq.=$_;
				  next ;
			   }elsif($list==1){
			      if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         print "$poss[$p]\n";
					 next A;
			      }
			   }elsif($fasta !=1){
				  print ;
			   }
			 }
		   }
		 }
		 close(DB);
	  }elsif($all==1){
		 print "\n";
		 open (DB, "$DB");
		 if(defined($SEQ_num_to_fetch)){ $num_sequence=$SEQ_num_to_fetch;
		 }else{ $num_sequence=@poss; }
		 A:for($p=0; $p < $num_sequence; $p++){
		   if($poss[$p]=~/\w*$target\w*/){
			 $matched=$possible{$poss[$p]};
			 seek(DB, ($matched-51), 0);
			 while(<DB>){
			   if($gdf==1){
			      if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         printf ("%-24s %-3d %-7d %-14s %4s\n", "$poss[$p]\/1\-$1", 1, $1, $poss[$p], '0.0');
					 next A;
			      }
			   }elsif(/^ {0,2}\/\// and $fasta==1){ # !!! DO NOT put $ in /^ {0,2}\/\// as there is something
				  $seq=~s/ //g;
				  if( ($SEQ_size_min < length($seq))&&(length($seq) < $SEQ_size_max) ){
					 print "\>$poss[$p]\n$seq\n"; $seq='';  next A;
				  }else{  $seq=''; $num_sequence++; next A; }
			   }elsif(/^ {0,2}\/\// and $fasta !=1){  # !!! DO NOT put $ in /^ {0,2}\/\// as there is something
				  print "\/\/\n";
				  next A;
			   }elsif(($fasta==1)&&(/^[\t ]+\w+/)){
				  $seq.=$_;
				  next ;
			   }elsif($list==1){
			      if(/ID[\t ]+$poss[$p] +\S+ +\S+ +(\d+)/){
			         printf "$poss[$p]\n";
					 next A;
			      }
			   }elsif($fasta !=1){
				  print ;
			   }
			 }
		   }
		 }
		 close(DB);
	  }

	  FASTA:
	  @poss = sort keys %possible;
	  if( (@poss >1)&&($all !=1)){
		 print "\n @poss","\n";
		 print chr(7);
		 print "\n There are more than a few seqs for $in[$i]";
		 print "\n be more specific! OR use -a option for all matched\n\n";
		 exit;
	  }elsif($all !=1){
		 print "\n";
		 open (FAS, "$FASTA");
		 B:for($p=0; $p < @poss; $p++){
		 if($poss[$p]=~/\w*$target\w*/){
			 $matched=$possible{$poss[$p]};
			 seek(FAS, ($matched-350), 0);
			 my $seq_found;
			 while(<FAS>){
			if((/^> {0,4}(\S+)/)&&($seq_found==1)){
				   next B;
	 			}elsif(/^> {0,4}($poss[$p])/){
				   print;
				   $seq_found=1;
				}elsif($seq_found==1){
				   print;
				}
			 }
		   }
		 }
		 close(FAS);
	  }elsif($all==1){
		 print "\n";
		 open (FAS, "$FASTA");
		 B2:for($p=0; $p < @poss; $p++){
		   if($poss[$p]=~/\w*$target\w*/){
			 $matched=$possible{$poss[$p]};
			 seek(FAS, ($matched-350), 0);
			 my $seq_found;
			 while(<FAS>){
				if((/^> {0,4}(\S+)/)&&($seq_found==1)){
				   next B2;
				}elsif(/^> *($poss[$p])/){
				   print;
				   $seq_found=1;
				}elsif($seq_found==1){
				   print;
				}
			 }
		   }
		 }
		 close(FAS);
	  }
   }
}



#__________________________________________________________________
# Title     : HELP_fetch_seq
# Usage     :
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Class     :
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
# Package   : Bio
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Scientist
# Version   : 1.0
# Used in   :
# Enclosed  :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub HELP_fetch_seq{
	 print "\n Usage: $0 [-cafghl] <any swissprot name entry>\n\n";
	 print "  -f is for FASTA output format only\n";
	 print "  -a is for ALL matched seq. \(same as using glob=\> *YEAST\)\n";
	 print "  -c is for Creating seq.idx file\n";
	 print "  -h is for HELP!\n";
	 print "  -g is for GDF file format output\n";
	 print "  -l is for list of match entries\(in 1 column\)\n";
	 print "  -s is for species option (input name mst be species (YEAST, RAT, HUMAN..)\n";
	 print "  n= is for Number of seq you want to get from swissprot\n";
	 print "  s= is for Size limit. Min seq size in swiss, s=10  -> minimum 11 aa seq.\n";
	 print "  S= is for Size limit. Max seq size in swiss, s=1000 -> get less than 1000 aa seq\n";
	 print "\n\n $0 uses ENV vars if exist for seq.dat and seq.idx";
	 print "\n     FETCHSWISS, SWDIR or SWISS to seq.dat or its path";
	 print "\n     SWINDEX or FETCHSWISSINDEX to seq.idx or its path\n";
	 print "\n   Option: set ENV vars, FASTADB, FASTAINDEX to your fasta db and fasta idx files\n";
	 print "\n   e.g. favourite fasta db and index file: seq.fa, seq.fa.idx in /usr/me/db/seq.fa\n";
	 print "\n        Your ENV set: setenv FASTADB=/usr/me/db/seq.fa\n";
	 print "\n        Your ENV set: setenv FASTAINDEX=/usr/me/db/seq.fa.idx\n";
	 print "\n   fetch also looks at your home dir and pwd for seq.dat, seq.idx \n";
	 print "\n   In  csh: \'setenv SWISS /your/DB/swiss/\'";
	 print "\n   In bash: \'export SWISS=/your/DB/swiss/\'\n";
	 print "   You can put absolute path for seq.dat and(or) seq.idx at prompt\n\n\n";
	 print chr(7);
}



#________________________________________________________________________
# Title     : read_dir_names_only
# Usage     : @all_dirs_list = @{&read_dir_names_only(\$absolute_path_dir_name, ....)};
# Function  : read any dir names and and then put in array.
# Example   :
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Todo      :
# Author    : A Biomatic
# Version   : 3.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_dir_names_only{
  my($in_dir, $i,$k, @possible_dirs,
	  @final_files, $full_dir, $pwd, $path,@read_files);
  $pwd=`pwd`; chomp($pwd); $full_dir=1;
  for($k=0; $k < @_; $k++){
	 if   ( ($_[$k] eq '.') || !(defined($_[$k]))){  $in_dir=$pwd;  }
	 elsif(!(ref($_[$k]))){   $in_dir=$_[$k];   }
	 elsif(ref($_[$k])){      $in_dir =${$_[$k]};    }
	 if($in_dir =~ /^([\w\-\.]+)$/){  $in_dir="$pwd\/$in_dir"; $full_dir = 0; }
	 else{ $full_dir =1; }
	 ##########  Main READING PART ##########
	 opendir(DIR1,"$in_dir");
	 @read_files = readdir(DIR1);
	 for($i=0; $i < @read_files; $i ++){
		$read_files[$i]="$in_dir\/$read_files[$i]";
		if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
		  $read_files[$i]=~s/\.\///; ## removing ./ in front of dirs (in bash)
		  push(@final_files, "$read_files[$i]");
		}
	 }
  }
  return([sort @final_files]);
}
