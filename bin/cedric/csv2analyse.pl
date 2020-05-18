#!/usr/bin/env perl
use Env;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use strict;
my @fam;
my @bucket;
my @aligner;
my @tree;
my %ignore;
my %full;
my %NN;
my %usedF;
my $N;
my %scoreH;
my %metricsH;


#$ignore{"parttreednd2"}=1;
#$ignore{"parttreednd2size"}=1;
#$ignore{"MAFFT-GINSI"}=1;
#$ignore{"MAFFT-FFTNS"}=1;
#$ignore{"ClustalO"}=1;

my $minseq=0;
my $maxseq=100000;
my $mrdelta=0.2;
my $score="sp";
my $metrics="ngap";
my $norm="original";
my $dir=getcwd;
my $reverse="no";


# Parse arguments
for (my $a=0; $a<=$#ARGV; $a++)
  {
    if ($ARGV[$a] eq "-minseq"){$minseq=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-maxseq"){$maxseq=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-mrdelta"){$mrdelta=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-score"){$score=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-metrics"){$metrics=$ARGV[++$a];}
	elsif ($ARGV[$a] eq "-norm"){$norm=$ARGV[++$a];}  # See normalization section
    elsif ($ARGV[$a] eq "-dir"){$dir=$ARGV[++$a];}
	elsif ($ARGV[$a] eq "-reverse"){$reverse=$ARGV[++$a];}
    else 
      {
	print "Unknown Flag: [$ARGV[$a]\n";
	die;
      }
  }


# Files
my %metricsF=csv2h("$dir/$metrics.csv");
my %scoreH=csv2h("$dir/$score.csv");
my %len=csv2h("$dir/len.csv");
#my %sp=csv2h("$dir/sp.csv");
#my %tc=csv2h("$dir/tc.csv");
#my %col=csv2h("$dir/col.csv");




# Info
$ignore{"Family"}=1;
$ignore{"nseq"}=1;
#my @fam=keys(%metricsF);
@fam=shrinklist(@fam);
@bucket=shrinklist(@bucket);
@aligner=shrinklist (@aligner);
@tree=shrinklist(@tree);


#@tree=("codnd", "parttreednd0");


# Normalization
foreach my $b (@bucket)
{
foreach my $aln (@aligner)
  {
    foreach my $f (@fam) 
      {
	foreach my $t1 (@tree)
	  {
		if ($norm ne "original")
		  {
			  if ($norm eq "normPerLen")  # metric / length
			    {
					my $a=$len{$f}{$b}{$aln}{$t1};
					if ($len{$f}{$b}{$aln}{$t1}>0  && $metricsF {$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}/$len{$f}{$b}{$aln}{$t1};}
					else {$metricsF{$f}{$b}{$aln}{$t1}="NA";}
			    }
			  elsif ($norm eq "normByLen") # metric * length
			    {
					if ($len{$f}{$b}{$aln}{$t1}>0  && $metricsF {$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}*$len{$f}{$b}{$aln}{$t1};}
					else {$metricsF{$f}{$b}{$aln}{$t1}="NA";}
				}
			  elsif ($norm eq "normPerSeq") # metric / nseq
			    {
					if ($len{$f}{nseq}{nseq}{nseq}>0 && $metricsF{$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}/$len{$f}{nseq}{nseq}{nseq};}
					else {$metricsF{$f}{$aln}{$t1}="NA";}
				}
			  elsif ($norm eq "normBySeq") # metric * nseq
			    {
					if ($len{$f}{nseq}{nseq}{nseq}>0 && $metricsF{$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}*$len{$f}{nseq}{nseq}{nseq};}
					else {$metricsF{$f}{$b}{$aln}{$t1}="NA";}
				}
			  elsif ($norm eq "normPerLenSeq")  # metric / length / nseq
			    {
					if ($len{$f}{nseq}{nseq}{nseq}>0 && $len{$f}{$b}{$aln}{$t1}>0  && $metricsF{$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}/$len{$f}{$b}{$aln}{$t1}/$len{$f}{nseq}{nseq}{nseq};}
					else {$metricsF{$f}{$b}{$aln}{$t1}="NA";}
				}
			  elsif ($norm eq "normByLenSeq")   # metric * length * nseq
			    {
					if ($len{$f}{nseq}{nseq}{nseq}>0 && $len{$f}{$b}{$aln}{$t1}>0  && $metricsF{$f}{$b}{$aln}{$t1} ne "NA"){$metricsH{$f}{$b}{$aln}{$t1}=$metricsF{$f}{$b}{$aln}{$t1}*$len{$f}{$b}{$aln}{$t1}*$len{$f}{nseq}{nseq}{nseq};}
					else {$metricsF{$f}{$b}{$aln}{$t1}="NA";}
				}
		  }
		else {%metricsH=%metricsF;}
	  }
      }
  }
}


my ($Tot, $TotN, $TotP);
my %fraction;
my @nfam;
my %pick;
foreach my $b (@bucket)
{
  my($BTot, $BTotN, $BTotP);
foreach my $aln (@aligner)
  {
    my ($ATot, $ATotN, $ATotP);
    foreach my $f (@fam) 
      {

	my ($tot, $totN, $totP);
	foreach my $t1 (@tree)
	  {
	    $full{$N}{metrics}=$metricsH{$f}{$b}{$aln}{$t1};
	    $full{$N}{score}=$scoreH{$f}{$b}{$aln}{$t1};
	    $full{$N}{aligner}=$aln;
	    
	    $N++;
	    
	    foreach my $t2 (@tree)
	      {
		
		if ($t1 eq $t2){next;}

		# Compute deltas
		my $use;
		my $h1=$metricsH{$f}{$b}{$aln}{$t1};
		my $h2=$metricsH{$f}{$b}{$aln}{$t2};
		my $D1=$h1-$h2;
		
		my $s1=$scoreH{$f}{$b}{$aln}{$t1};
		my $s2=$scoreH{$f}{$b}{$aln}{$t2};
		my $D2=$s1-$s2;
		
		
		$NN{$b}{$aln}{$f}{tot}++;
		if ( $s1 eq "NA" || $s2 eq "NA" || $h1 eq "NA" || $h2 eq "NA")
		  {
		    $NN{$b}{$aln}{$f}{unused}++;
		  }
		elsif ($len{$f}{nseq}{nseq}{nseq}<$minseq || $len{$f}{nseq}{nseq}{nseq}>$maxseq)
		  {
		    $NN{$b}{$aln}{$f}{unused}++;
		  }
		elsif (($h1==$h2 && $h1==0)||(abs($D1)*2)/($h1+$h2)<$mrdelta)
		  {
		    $NN{$b}{$aln}{$f}{unused}++;
		  }
		elsif (abs($D1)<=0.001 && abs($D2)<=0.001)
		  {
		    $NN{$b}{$aln}{$f}{unused}++;
		  }
		# elsif ($norm ne "original" && abs($D1)<=0.0000000000001 && abs($D2)<=0.0000000000001)
		# #elsif ($norm ne "original" && abs($D1)<=0.001 && abs($D2)<=0.001)
		#   {
		#     $NN{$b}{$aln}{$f}{unused}++;
		#   }
		elsif (($D2>=0 && $D1>=0) || ($D2<=0 && $D1<=0))
		  {
			#printf "+usedfamily $aln $f\n";
		    $usedF{$b}{$aln}{$f}=1;  # Used family
			$usedF{$b}{all}{$f}=1;
			$usedF{all}{all}{$f}=1;

		    $Tot++;
		    $TotP++;

			$BTot++;
			$BTotP++;
		    
		    $ATot++;
		    $ATotP++;
		    
		    $tot++;
		    $totP++;

		    $use=1;
		  }
	 	elsif (($D2>0 && $D1<0) || ($D2<0 && $D1>0))
		  {
		    $usedF{$b}{$aln}{$f}=1;  # Used family
			$usedF{$b}{all}{$f}=1;
			$usedF{all}{all}{$f}=1;

		    $Tot++;
		    $TotN++;
		    
			$BTot++;
			$BTotN++;

		    $ATot++;
		    $ATotN++;
		    
		    $tot++;
		    $totN++;
		    $use=1;
		  }
		if ($use) 
		  {
		    if ($h1>$h2)
		      {
				if ($reverse eq "no")
				{
				$pick{$b}{$aln}{maxS}+=$s2;
				$pick{$b}{$aln}{maxN}++;
				$pick{$b}{$aln}{minS}+=$s1;
				$pick{$b}{$aln}{minN}++;
				}
				else
				{
				$pick{$b}{$aln}{maxS}+=$s1;
				$pick{$b}{$aln}{maxN}++;
				$pick{$b}{$aln}{minS}+=$s2;
				$pick{$b}{$aln}{minN}++;
				}
			}
		    else
		      {
				if ($reverse eq "no")
				{
				$pick{$b}{$aln}{maxS}+=$s1;
				$pick{$b}{$aln}{maxN}++;
				$pick{$b}{$aln}{minS}+=$s2;
				$pick{$b}{$aln}{minN}++;
				}
				else
				{
				$pick{$b}{$aln}{maxS}+=$s2;
				$pick{$b}{$aln}{maxN}++;
				$pick{$b}{$aln}{minS}+=$s1;
				$pick{$b}{$aln}{minN}++;
				}
		      }
		  }
			
		
	      }  # Finish t2
	  }      # Finish t1

	$NN{$b}{$aln}{all}{tot}+=$NN{$b}{$aln}{$f}{tot};
	$NN{$b}{$aln}{all}{unused}+=$NN{$b}{$aln}{$f}{unused};

	if ($reverse eq "no")
		{
		my $R1=($tot>0)?$totN/$tot:0;	# Number of points in N / total used point
		my $R2=($NN{$b}{$aln}{$f}{tot}>0)?$NN{$b}{$aln}{$f}{unused}/$NN{$b}{$aln}{$f}{tot}:0;   # Number of unused / total trees**2
		#printf ("FAM::$f: [$b] [$aln]: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d)\n", $R1,$totN,$tot, $R2,$NN{$b}{$aln}{$f}{unused},$NN{$b}{$aln}{$f}{tot});  
		printf ("[FAM:BUCKET:ALIGNER]:[$f:$b:$aln]:: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d)\n", $R1,$totN,$tot, $R2,$NN{$b}{$aln}{$f}{unused},$NN{$b}{$aln}{$f}{tot});
		}
	else
		{
		my $R1=($tot>0)?$totP/$tot:0;
		my $R2=($NN{$aln}{$f}{tot}>0)?$NN{$aln}{$f}{unused}/$NN{$aln}{$f}{tot}:0;   # Number of unused / total trees**2
		#printf ("FAM::$f: [$b] [$aln]: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d)\n", $R1,$totP,$tot, $R2,$NN{$b}{$aln}{$f}{unused},$NN{$b}{$aln}{$f}{tot});  
		printf ("[FAM:BUCKET:ALIGNER]:[$f:$b:$aln]:: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d)\n", $R1,$totP,$tot, $R2,$NN{$b}{$aln}{$f}{unused},$NN{$b}{$aln}{$f}{tot});
		}  	 
      }  # Finish f
    
    $NN{$b}{all}{all}{tot}+=$NN{$b}{$aln}{all}{tot};
    $NN{$b}{all}{all}{unused}+=$NN{$b}{$aln}{all}{unused};

	$pick{$b}{all}{maxS}+=$pick{$b}{$aln}{maxS};
	$pick{$b}{all}{maxN}+=$pick{$b}{$aln}{maxN};
	$pick{$b}{all}{minS}+=$pick{$b}{$aln}{minS};
	$pick{$b}{all}{minN}+=$pick{$b}{$aln}{minN};
    

	if ($reverse eq "no")
		{
		my $R1=($ATot>0)?$ATotN/$ATot:0;
		my $R2=($NN{$b}{$aln}{all}{tot}>0)?$NN{$b}{$aln}{all}{unused}/$NN{$b}{$aln}{all}{tot}:0;
		my @used_families=keys (%{$usedF{$b}{$aln}});
		my $R3=$#used_families+1;
		my $min_acc=($pick{$b}{$aln}{minN}>0)?$pick{$b}{$aln}{minS}/$pick{$b}{$aln}{minN}:0;  # TC or SP score on average if the worst MSA is picked every time
		my $max_acc=($pick{$b}{$aln}{maxN}>0)?$pick{$b}{$aln}{maxS}/$pick{$b}{$aln}{maxN}:0;  # TC or SP score on average if the best MSA is picked every time
		#printf ("ALN::$aln [$b] : %.3f (%d/%d) -- UNUSED: %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1, $ATotN,$ATot, $R2, $NN{$b}{$aln}{all}{unused},$NN{$b}{$aln}{all}{tot}, $min_acc, $max_acc, $max_acc-$min_acc);
		printf ("[BUCKET:ALIGNER]:[$b:$aln]:: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1, $ATotN,$ATot, $R2, $NN{$b}{$aln}{all}{unused},$NN{$b}{$aln}{all}{tot}, $min_acc, $max_acc, $max_acc-$min_acc);
		}
		else
		{
		my $R1=($ATot>0)?$ATotP/$ATot:0;
		my $R2=($NN{$b}{$aln}{all}{tot}>0)?$NN{$b}{$aln}{all}{unused}/$NN{$b}{$aln}{all}{tot}:0;
		my @used_families=keys (%{$usedF{$b}{$aln}});
		my $R3=$#used_families+1;
		my $min_acc=($pick{$b}{$aln}{minN}>0)?$pick{$b}{$aln}{minS}/$pick{$b}{$aln}{minN}:0;  # TC or SP score on average if the worst MSA is picked every time
		my $max_acc=($pick{$b}{$aln}{maxN}>0)?$pick{$b}{$aln}{maxS}/$pick{$b}{$aln}{maxN}:0;  # TC or SP score on average if the best MSA is picked every time
		#printf ("ALN::$aln [$b] : %.3f (%d/%d) -- UNUSED: %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1, $ATotP,$ATot, $R2, $NN{$b}{$aln}{all}{unused},$NN{$b}{$aln}{all}{tot}, $min_acc, $max_acc, $max_acc-$min_acc);
		printf ("[BUCKET:ALIGNER]:[$b:$aln]:: %.3f (%d/%d) -- UNUSED: %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1, $ATotP,$ATot, $R2, $NN{$b}{$aln}{all}{unused},$NN{$b}{$aln}{all}{tot}, $min_acc, $max_acc, $max_acc-$min_acc);
		}
    
  }  # finish aln

    $NN{all}{all}{all}{tot}+=$NN{$b}{all}{all}{tot};
    $NN{all}{all}{all}{unused}+=$NN{$b}{all}{all}{unused};

	$pick{all}{all}{maxS}+=$pick{$b}{all}{maxS};
	$pick{all}{all}{maxN}+=$pick{$b}{all}{maxN};
	$pick{all}{all}{minS}+=$pick{$b}{all}{minS};
	$pick{all}{all}{minN}+=$pick{$b}{all}{minN};

	if ($reverse eq "no")
	{
	my $R1=($BTot>0)?$BTotN/$BTot:0;
	my $R2=($NN{$b}{all}{all}{tot}>0)?$NN{$b}{all}{all}{unused}/$NN{$b}{all}{all}{tot}:0;
	my @used_families=keys (%{$usedF{$b}{all}});
	my $R3=$#used_families+1;
	my $min_acc=($pick{$b}{all}{minN}>0)?$pick{$b}{all}{minS}/$pick{$b}{all}{minN}:0;
	my $max_acc=($pick{$b}{all}{maxN}>0)?$pick{$b}{all}{maxS}/$pick{$b}{all}{maxN}:0;
	#printf ("GLOBAL[$b]: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$BTotN,$BTot, $R2,$NN{$b}{all}{all}{unused},$NN{$b}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
	printf ("[BUCKET]:[$b]:: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$BTotN,$BTot, $R2,$NN{$b}{all}{all}{unused},$NN{$b}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
	}
	else
	{
	my $R1=($BTot>0)?$BTotP/$BTot:0;
	my $R2=($NN{$b}{all}{all}{tot}>0)?$NN{$b}{all}{all}{unused}/$NN{$b}{all}{all}{tot}:0;
	my @used_families=keys (%{$usedF{$b}{all}});
	my $R3=$#used_families+1;
	my $min_acc=($pick{$b}{all}{minN}>0)?$pick{$b}{all}{minS}/$pick{$b}{all}{minN}:0;
	my $max_acc=($pick{$b}{all}{maxN}>0)?$pick{$b}{all}{maxS}/$pick{$b}{all}{maxN}:0;
	#printf ("GLOBAL[$b]: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$BTotP,$BTot, $R2,$NN{$b}{all}{all}{unused},$NN{$b}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
	printf ("[BUCKET]:[$b]:: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$BTotP,$BTot, $R2,$NN{$b}{all}{all}{unused},$NN{$b}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
	}
}  # Finish b


if ($reverse eq "no")
{
my $R1=($Tot>0)?$TotN/$Tot:0;
my $R2=($NN{all}{all}{all}{tot}>0)?$NN{all}{all}{all}{unused}/$NN{all}{all}{all}{tot}:0;
my @used_families=keys (%{$usedF{all}{all}});
my $R3=$#used_families+1;
my $min_acc=($pick{all}{all}{minN}>0)?$pick{all}{all}{minS}/$pick{all}{all}{minN}:0;
my $max_acc=($pick{all}{all}{maxN}>0)?$pick{all}{all}{maxS}/$pick{all}{all}{maxN}:0;
printf ("ALL: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$TotN,$Tot, $R2,$NN{all}{all}{all}{unused},$NN{all}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
}
else
{
my $R1=($Tot>0)?$TotP/$Tot:0;
my $R2=($NN{all}{all}{all}{tot}>0)?$NN{all}{all}{all}{unused}/$NN{all}{all}{all}{tot}:0;
my @used_families=keys (%{$usedF{all}{all}});
my $R3=$#used_families+1;
my $min_acc=($pick{all}{all}{minN}>0)?$pick{all}{all}{minS}/$pick{all}{all}{minN}:0;
my $max_acc=($pick{all}{all}{maxN}>0)?$pick{all}{all}{maxS}/$pick{all}{all}{maxN}:0;
printf ("ALL: %.3f (%d/%d) -- UNUSED %.3f (%d/%d) [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1,$TotP,$Tot, $R2,$NN{all}{all}{all}{unused},$NN{all}{all}{all}{tot},$min_acc, $max_acc, $max_acc-$min_acc);
}



exit (0);


  

sub csv2h
	 {
	   #use Tie::IxHash;
	   #tie my %h, 'Tie::IxHash';	
	   my %h;
	   my ($f)=@_;   
	   my $n;
	   if ( !-e $f || !($f=~/.*\.csv/)){return %h;}
	   
	   print "----Process $f\n";
	   open (F, "$f");
	   while (<F>)
	     {
	       my $l=$_;
	       my @ll;
	       my @lll=split (/\,/,$l);
		   
	      
	       foreach my $x (@lll)
		 {
		   $x=~s/\r//g;
		   $x=~s/\n//g;
		   $x=~s/\"//g;
		   @ll=(@ll,$x);
		 }
	       if ($n>2)  # Values
		 {
		   push(@fam,$ll[0]);

		   for (my $a=1;$a<=$#aligner; $a++)  # For each column, "family" excluded
		     {
			   $ll[$a]=$lll[$a];
		       #$ll[$a]=~s/\"//g;   # Value at column $a    # Check this line;; incorrect value
		       $h{$ll[0]}{$bucket[$a]}{$aligner[$a]}{$tree[$a]}=$ll[$a];    # h{family}{bucket}{aligner}{tree}=value
		       #if ($ll[$a] eq "NA"){$ignore{$ll[0]}=1;print "***** $ll[0]\n";}

		     }
		 }
	       elsif ($n==0)
		 {
		   @bucket=@ll;   # Bucket header
		 }
	       elsif ($n==1)
		 {
		   @aligner=@ll;  # Aligner header
		 }
		   elsif ($n==2)
		 {
		   @tree=@ll;     # Tree header
		 }
	       $n++;
	     }
	   close (F);
	   return %h;

	 }


	 


sub shrinklist 
	   {
	   use Tie::IxHash;
	   my @l=@_;
	   tie my %h, 'Tie::IxHash';

	   foreach my $e (@l){if ($e ne "" && $e ne "Family" && !$ignore{$e}){$h{$e}=1;}}
	   #return sort keys (%h);
	   return keys (%h);
	 }
	   
