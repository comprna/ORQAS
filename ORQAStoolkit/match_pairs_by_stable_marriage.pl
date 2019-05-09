#!/usr/bin/perl -w

use strict;

######################################################################
# Written by Eduardo Eyras
# You may distribute this module under the same terms as perl itself
# Let the code begin...

my $verbose = 0;

my ($file) = @ARGV;
unless ($file){
    print STDERR "Use: $0 file\n";
    exit(0);
}

my %score;                  # score for every pair of proteins
my %potential_partners1;    # potential partners for proteins in list 1
my %potential_partners2;    # potential partners for proteins in list 2

my %prots1; # proteins in gene 1
my %prots2; # proteins in gene 2

my %married_prot1;
my %married_prot2;

my %partner_in_1;
my %partner_in_2;

my %gene_partner;

my %seq;

open(IN,"<$file") or die("cannot open $file as input");
while(<IN>){
    chomp;
    # ENSG00000175806     ENSMUSG00000054733  0.46153846153846156 MAVFGMGCFWGAER... MAVFGMGCFWGAER...

    my ($gene1, $gene2, $score, $seq1, $seq2) = split;

    my $prot1 = $seq1;
    my $prot2 = $seq2;

    $seq{$prot1} = $seq1 if $seq1;
    $seq{$prot2} = $seq2 if $seq2;

    # for each gene pair, we store these pair of proteins:

    $gene_partner{$gene1} = $gene2;
    
    $prots1{$gene1}{$prot1}++;
    $prots2{$gene2}{$prot2}++;
    push( @{$potential_partners1{$prot1}}, $prot2 );    
    push( @{$potential_partners2{$prot2}}, $prot1 );    

    $score{$prot1}{$prot2} = $score;
}

############################################################
# this method takes all the pairs with the corresponding scores
# and finds the best pairs using the 'stable-marriage' algorithm.

# The Algorithm has the following condition for solution:
# there are no two elements that they're not paired-up to each other but
# they have better score with each other ( $score_matrix is higher ) 
# than with their current corresponding partners

# the main different between this optimization algorithm
# and a 'best-reciprocal-pairs' approach is that
# 'stable-marriage' produces rather 'best-available-pairs', so it keeps matches which
# not being maximal are still optimal. It warranties only on pair
# per element, and this is the best available one, i.e. ' you like
# C. Schiffer but she likes D. Copperfield more than she likes you so
# you have to stick to the next one in your priority list if available'.


# we iterate over all the genes in species 1:
foreach my $gene1 (keys %prots1 ){
    
    # we iterate over all proteins in gene1
    my @prots1 = keys %{$prots1{$gene1}};
    my @unmarried_ones = @prots1;
    
 MARRIAGE:
  while ( @unmarried_ones ){
      
      # pick one of them
      my $prot1 = shift @unmarried_ones;
    
      # get the potential partners
      my @partners1 = @{$potential_partners1{$prot1}};
    
      # sort them in descending order according to score:
      my @sorted_partners_in_2 = sort_scores( $prot1, \@partners1 );
  
      # go over the partners until you get married or run out of partners  
    PARTNER:
      while( @sorted_partners_in_2 && !defined($married_prot1{$prot1}) ){
	  
	  print "checking partner list for $prot1\n" if $verbose;
	  my $potential_partner_in_list2 = shift( @sorted_partners_in_2 );
	  
	  print "looking at $potential_partner_in_list2\n" if $verbose;
	  # check whether it is already married
	  if ( $married_prot2{ $potential_partner_in_list2 } 
	       && $married_prot2{ $potential_partner_in_list2 } == 1 ){
	      
	      # is it married to another target?
	      if ( $partner_in_1{$potential_partner_in_list2} 
		   &&  !( $partner_in_1{ $potential_partner_in_list2 } eq $prot1 ) ){
		  
		  # is it a 'worse' marriage?
		  if ( $score{$partner_in_1{$potential_partner_in_list2}}{$potential_partner_in_list2}
		       < $score{$prot1}{$potential_partner_in_list2} ){
		      
		      # put the divorced one back into the pool only if it has more potential partners
		      if ( @{$potential_partners1{$partner_in_1{$potential_partner_in_list2}}} ){
			  push ( @unmarried_ones, $partner_in_1{$potential_partner_in_list2} );
		      }
		      
		      # divorce the 'worse partner'
		      print "divorcing ".$partner_in_1{$potential_partner_in_list2}."\n" if $verbose;
		      delete $married_prot1{ $partner_in_1{ $potential_partner_in_list2 } };
		      delete $partner_in_2{ $partner_in_1{ $potential_partner_in_list2 } };
		      delete $partner_in_1{ $potential_partner_in_list2 };
		      
		      # let be happier marriage
		      $married_prot1{ $prot1 } = 1;
		      $married_prot2{ $potential_partner_in_list2 } = 1;
		      $partner_in_1{ $potential_partner_in_list2 } = $prot1;
		      $partner_in_2{ $prot1 } = $potential_partner_in_list2;
		      print "new marriage: $prot1 -- $potential_partner_in_list2\n" if $verbose;
		      next MARRIAGE;
		      
		  }
		  else{
		      # look at the next potential partner in list2
		      next PARTNER;
		  }
	      }
	      # hmm, this prot2 (in list2)  is married, to whom?
	      elsif ( $partner_in_1{ $potential_partner_in_list2 } eq $prot1 ) {
		  # hey, we have already a happy couple
		  $partner_in_2{ $prot1 } = $potential_partner_in_list2;
		  next MARRIAGE;
	      }
	      elsif ( !defined( $partner_in_1{ $potential_partner_in_list2 } ) ){
		  # we have a cheater!
		  $married_prot2{ $potential_partner_in_list2 } = 0;
		  next PARTNER;
	      }
	  }
	  else{
	      
	      # this prot2 ( in list 2 ) is still single, let be marriage:
	      $married_prot1{ $prot1 } = 1;
	      $married_prot2{ $potential_partner_in_list2 } = 1;
	      $partner_in_1{ $potential_partner_in_list2 } = $prot1;
	      $partner_in_2{ $prot1 } = $potential_partner_in_list2;
	      print "setting partner{ $prot1 } = $potential_partner_in_list2\n" if $verbose;
	      next MARRIAGE;
	  }
      } # end of PARTNER
  }   # end of MARRIAGE
    
    foreach my $prot1 (@prots1){
	if ( $partner_in_2{$prot1} ){
	    my $seq1 = "";
	    my $seq2 = "";
	    $seq1 = $seq{$prot1} if $seq{$prot1};
	    $seq2 = $seq{$partner_in_2{$prot1}} if $seq{$partner_in_2{$prot1}};
	    my $s = join "\t", ($gene1, $gene_partner{$gene1}, 
				$score{$prot1}{$partner_in_2{$prot1}}, $seq1, $seq2);
	    print $s."\n";
	}
    }
}


##################################
#





############################################################
# sort the potential partners by score in descending order
sub sort_scores{
    my ($prot1, $partners) = @_;
    my @partners = @$partners;
    my @sorted_partners =  map {$_->[1]} sort {$b->[0] <=> $a->[0]}  map { [$score{$prot1}{$_}, $_] } @partners;
    return @sorted_partners;
}

############################################################




