#!/usr/bin/env perl
use Env;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use strict;
my @tree;
my @aligner;
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
my $dir=getcwd;


for (my $a=0; $a<=$#ARGV; $a++)
  {
    if ($ARGV[$a] eq "-minseq"){$minseq=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-maxseq"){$maxseq=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-mrdelta"){$mrdelta=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-score"){$score=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-metrics"){$metrics=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-dir"){$dir=$ARGV[++$a];}
    else 
      {
	print "Unknown Flag: [$ARGV[$a]\n";
	die;
      }
  }

my %ngap=csv2h("$dir/ngap.csv");
my %ngap2=csv2h("$dir/ngap2.csv");
my %homo=csv2h("$dir/homo.csv");
my %whomo=csv2h("$dir/w_homo.csv");
my %whomo2=csv2h("$dir/w_homo2.csv");
my %Nngap;#defined below
my %whomo3;#defined below
my %whomo4;#defined below

my %len=csv2h("$dir/len.csv");
my %sp=csv2h("$dir/sp.csv");
my %tc=csv2h("$dir/tc.csv");
my %col=csv2h("$dir/col.csv");


$ignore{"Family"}=1;
$ignore{"nseq"}=1;
my @fam=keys(%len);
@fam=shrinklist(@fam);
@aligner=shrinklist (@aligner);
@tree=shrinklist(@tree);



#@tree=("codnd", "parttreednd0");


foreach my $aln (@aligner)
  {
    foreach my $f (@fam) 
      {
	foreach my $t1 (@tree)
	  {
	    if ($len{$f}{$aln}{$t1}>0  && $Nngap {$f}{$aln}{$t1} ne "NA"){$Nngap {$f}{$aln}{$t1}=$ngap  {$f}{$aln}{$t1}/$len{$f}{$aln}{$t1};}
	    if ($len{$f}{$aln}{$t1}>0  && $whomo2{$f}{$aln}{$t1} ne "NA"){$whomo3{$f}{$aln}{$t1}=$whomo2{$f}{$aln}{$t1}*$len{$f}{$aln}{$t1};}
	    if ($len{$f}{nseq}{nseq}>0 && $whomo2{$f}{$aln}{$t1} ne "NA"){$whomo4{$f}{$aln}{$t1}=$whomo2{$f}{$aln}{$t1}*$len{$f}{nseq}{nseq};}
	    
	    else
	      {
		$Nngap{$f}{$aln}{$t1}="NA";
		$whomo3{$f}{$aln}{$t1}="NA";
		$whomo4{$f}{nseq}{nseq}="NA";
	      }
	  }
      }
  }

#select score
if    ($score eq "sp") {%scoreH=%sp;}
elsif ($score eq "col"){%scoreH=%col;}
elsif ($score eq "tc") {%scoreH=%tc;}

#select metrics
if    ($metrics eq "ngap")  {%metricsH=%ngap;}
elsif ($metrics eq "ngap2") {%metricsH=%ngap2;}
elsif ($metrics eq "len")   {%metricsH=%len;}
elsif ($metrics eq "Nngap") {%metricsH=%Nngap;}

elsif ($metrics eq "homo")  {%metricsH=%homo;}
elsif ($metrics eq "whomo") {%metricsH=%whomo;}
elsif ($metrics eq "whomo2"){%metricsH=%whomo2;}
elsif ($metrics eq "whomo3"){%metricsH=%whomo3;}
elsif ($metrics eq "whomo4"){%metricsH=%whomo4;}




my ($Tot, $TotN, $TotP);
my %fraction;
my @nfam;
my %pick;
foreach my $aln (@aligner)
  {
    my ($ATot, $ATotN, $ATotP);
    foreach my $f (@fam) 
      {

	my ($tot, $totN, $totP);
	foreach my $t1 (@tree)
	  {
	    $full{$N}{metrics}=$metricsH{$f}{$aln}{$t1};
	    $full{$N}{score}=$scoreH{$f}{$aln}{$t1};
	    $full{$N}{aligner}=$aln;
	    
	    $N++;
	    
	    foreach my $t2 (@tree)
	      {

		# Compute deltas
		my $use;
		my $h1=$metricsH{$f}{$aln}{$t1};
		my $h2=$metricsH{$f}{$aln}{$t2};
		my $D1=$h1-$h2;
		
		my $s1=$scoreH{$f}{$aln}{$t1};
		my $s2=$scoreH{$f}{$aln}{$t2};
		my $D2=$s1-$s2;
		
		
		$NN{$aln}{$f}{tot}++;
		if ( $s1 eq "NA" || $s2 eq "NA" || $h1 eq "NA" || $h2 eq "NA")
		  {
		    $NN{$aln}{$f}{unused}++;
		  }
		elsif ($len{$f}{nseq}{nseq}<$minseq || $len{$f}{nseq}{nseq}>$maxseq)
		  {
		    $NN{$aln}{$f}{unused}++;
		  }
		elsif (($h1==$h2 && $h1==0)||(abs($D1)*2)/($h1+$h2)<$mrdelta)
		  {
		    $NN{$aln}{$f}{unused}++;
		  }
		elsif ( abs($D1)<=0.001 && abs($D2)<=0.001)
		  {
		    $NN{$aln}{$f}{unused}++;
		  }
		elsif (($D2>0 && $D1>=0) || ($D2<0 && $D1<=0))
		  {
		    $usedF{$aln}{$f}=1;
		    $Tot++;
		    $TotP++;
		    
		    $ATot++;
		    $ATotP++;
		    
		    $tot++;
		    $totP++;

		    $use=1;
		  }
	 	elsif (($D2>0 && $D1<=0) || ($D2<0 && $D1>=0))
		  {
		    $usedF{$aln}{$f}=1;
		    $Tot++;
		    $TotN++;
		    
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
			$pick{$aln}{maxS}+=$s2;
			$pick{$aln}{maxN}++;
			$pick{$aln}{minS}+=$s1;
			$pick{$aln}{minN}++;
			}
		    else
		      {
			$pick{$aln}{maxS}+=$s1;
			$pick{$aln}{maxN}++;
			$pick{$aln}{minS}+=$s2;
			$pick{$aln}{minN}++;
		      }
		  }
			
		
	      }
	  }
	$NN{$aln}{all}{tot}+=$NN{$aln}{$f}{tot};
	$NN{$aln}{all}{unused}+=$NN{$aln}{$f}{unused};
	
	my $R1=($tot>0)?$totN/$tot:0;   # Number of points in N / total used
	my $R2=($NN{$aln}{$f}{tot}>0)?$NN{$aln}{$f}{unused}/$NN{$aln}{$f}{tot}:0;   # Number of unused / total
	printf ("FAM::$f: [$aln]: %.3f -- UNUSED: %.3f\n", $R1, $R2);   
      }
    
    $NN{all}{all}{tot}+=$NN{$aln}{all}{tot};
    $NN{all}{all}{unused}+=$NN{$aln}{all}{unused};
    
    my $R1=($ATot>0)?$ATotN/$ATot:0;
    my $R2=($NN{$aln}{all}{tot}>0)?$NN{$aln}{all}{unused}/$NN{$aln}{all}{tot}:0;
    my @used_families=keys (%{$usedF{$aln}});
    my $R3=$#used_families+1;
    my $min_acc=$pick{$aln}{minS}/$pick{$aln}{minN};  # TC or SP score on average if the worst MSA is picked every time
    my $max_acc=$pick{$aln}{maxS}/$pick{$aln}{maxN};  # TC or SP score on average if the best MSA is picked every time
    
    printf ("ALN::$aln : %.3f -- UNUSED: %.3f [$R3 Families][Min Acc %.3f Max Acc %.3f Delta Acc %.3f]\n", $R1, $R2, $min_acc, $max_acc, $max_acc-$min_acc);
  }
my $R1=$TotN/$Tot;
my $R2=($NN{all}{all}{tot}>0)?$NN{all}{all}{unused}/$NN{all}{all}{tot}:0;
printf ("GLOBAL: %.3f -- UNUSED %.3f\n", $R1, $R2);

exit (0);


  

sub csv2h
	 {
	   my ($f)=@_;
	   my %h;	   
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
	       
	       if ($n>1)
		 {
		   for (my $a=1;$a<=$#aligner; $a++)
		     {
		       $ll[$a]=~s/\"//g;
		       $h{$ll[0]}{$aligner[$a]}{$tree[$a]}=$ll[$a];
		       #if ($ll[$a] eq "NA"){$ignore{$ll[0]}=1;print "***** $ll[0]\n";}
		       
		     }
		 }
	       elsif ($n==0)
		 {
		   @aligner=@ll;
		 }
	       elsif ($n==1)
		 {
		   @tree=@ll;
		 }
	       $n++;
	     }
	   close (F);
	   return %h;

	 }


	 


sub shrinklist 
	   {
	   my @l=@_;
	   my %h;

	   foreach my $e (@l){if ($e ne "" && $e ne "Family" && !$ignore{$e}){$h{$e}=1;}}
	   return sort keys (%h);
	 }
	   
