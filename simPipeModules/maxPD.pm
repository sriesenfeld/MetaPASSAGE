#!/usr/bin/perl
## S.J. Riesenfeld, Dongying Wu
## Updated: Oct 2010
## MetaPASSAGE module: maxPD
##
## Module created by Riesenfeld (with only very minor changes)
##   from Dongying Wu's maxPD.pl script
## 

=head2 maxPD_list
 
 Usage: maxPD_list(%args)
 Function: Given a phylogenetic tree, computes a list of taxa forming
   a subtree with maximum phylogenetic diversity (PD)
 Note: Taxa cannot have :,;()'"| or space in the name.
 Example: 
   # Assume $tree is a file containing a phylogenetic tree with
   # taxonomic information in leaf labels
   maxPD_list('-t', $tree, '-n' 10, '-o', 'output');
 Returns: A file giving a list of taxa (with corresponding information
   about contributions to PD)
 Args: A hash %args where the keys and values are as described
   below. The keys '-t', '-n', and '-o' must be defined; the rest are
   optional.
   Key: '-t'; Value: File containing a phylogenetic tree.
   Key: '-n'; Value: Positive integer giving the number of taxa to be
     selected.
   Key: '-o'; Value: Name of output file.

   Key: '-l'; Value: Name of file containing a list of taxa that are
     required to belong in the output. The function will start from
     there and pick the other taxa.  Key: '-pool'; Value: Name of file
     containing a list of taxa that are the only taxa permitted to
     appear in the output; the script will not choose any taxa outside
     this list.

=cut

package maxPD;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(maxPD_list);

sub maxPD_list(@) {
    my %opt=@_;
    my $tree_file=$opt{'-t'};
    my $input_number=$opt{'-n'};
    my $output_list_file=$opt{'-o'};    
    if (!defined ($tree_file) or !defined($input_number) or !defined($output_list_file)) {
	die "Missing required options for maxPD_list! ('-t', '-n', and '-o' are required)\n";
    }
    my $input_list=$opt{'-l'};
    my $gml_out=$opt{'-gml'};
    my $in_tree;    
    my $gml_out_file=$output_list_file.".gml";
    my $pool=$opt{'-pool'};
    my %pool;
    my $first_pool;

    if($pool){
	open(POOL,$pool)||die "cannot open pool file $pool";
	while(<POOL>){
	    $_=~s/\|/__/g; ##take care of "|"
	    $_=~s/\"//g; ##take care of """
	    $_=~s/\'//g; ##take care of "'"
	    my ($acc)=split(/\s+/);
	    if(!$first_pool){$first_pool=$acc;} 
	    $pool{$acc}=1;
	}
	close POOL;
    }
    
    if(!($input_number >= 1)){
	die "Please give get_maxPD a number >=1 \n";
    }
    
    my @input_list;
    
    if($input_list){
	open (INLIST,$input_list) || die "cannot open input list $input_list\n";
	while(<INLIST>){
	    $_=~s/\|/__/g; ##take care of "|"
	    $_=~s/\"//g; ##take care of """
	    $_=~s/\'//g; ##take care of "'"
	    my @t=split(/\s+/);
	    push(@input_list,$t[0]);
	}
	close INLIST;
    }
    
    
    open(IN,$tree_file) || die "cannot open tree file $tree_file\n";
    while(<IN>){
	$in_tree.=$_;
    }
    close IN;
    open (OUT,">$output_list_file") || die "Please define output file\n"; 
    my $tree_obj=DWU_TREE->new($in_tree);
    if($pool){
	$tree_obj->pool_only(\%pool);
    }
    
    my $max_list;
    
    if(@input_list){
	foreach my $acc(@input_list){
	    print OUT "# ".$acc."\n";
	}
	$tree_obj->build_sub_tree(@input_list);
	$tree_obj->get_max_PD($input_number,*OUT);
    }
    else{
	my $start_acc;
	if($first_pool){
	    $start_acc=$first_pool;
	}
	else{
	    $start_acc=$tree_obj->initial_seed();
	}
	print OUT $start_acc."\n";
	$tree_obj->build_sub_tree($start_acc);
	$input_number--;
    if($input_number >= 1){
	$tree_obj->get_max_PD($input_number,*OUT);
    }
    }
    
    close OUT;
    
    
    if(defined($gml_out) and ($gml_out=~/^[1yY]/)){
	$tree_obj->output_gml($gml_out_file);
    }
}

1;

package DWU_TREE;

sub output_gml {
    my ($class,$output_gml)=@_;
    open (GML,">$output_gml") || die "cannot output to $output_gml\n";
    print GML "graph [ \n";
    my %gml_id;
    my $gml_id=0;
    foreach my $key (keys %{$class->{matrix}}){
	$gml_id++;
        $gml_id{$key}=$gml_id;

        print GML "node [id ".$gml_id." ";
        if (($key=~/^T\d+/)&&($class->{sub_tree}->{$key})){
           if($class->{sub_tree}->{$key} ==2){ 
       	    print GML "lable \"".$class->{id_tn}->{$key}."\" graphics [ w 50 h 50 type \"circle\" fill \"\#FF0000\"]]\n";
	}
	   else{
            print GML "lable \"".$class->{id_tn}->{$key}."\" graphics [ w 50 h 50 type \"circle\" fill \"\#0000FF\"]]\n";
	   }
	}
        elsif($key=~/^T\d+/){
            print GML "lable \""."\" graphics [ w 50 h 50 type \"circle\" fill \"\#666666\"]]\n";
	}
        else{
            print GML "lable \""."\" graphics [ w 3 h 3 type \"circle\" fill \"\#666666\"]]\n";
	}
    }

    my %inlink;
   
    foreach my $from(keys %{$class->{matrix}}){
    foreach my $to(keys %{$class->{matrix}->{$from}}) {
    
	my ($source,$target)=sort($from,$to);

          my $link=$source." ".$target;
          if($inlink{$link}){
	      next;
	  }
          $inlink{$link}=1;          
          my $line_color="666666";
          if($class->{sub_tree}->{$source} && $class->{sub_tree}->{$target}){
          $line_color="FF0000";      
	  }
          elsif($class->{matrix}->{$source}->{$target} <= 0.0001){
          $line_color="FFFFFF";   
          } 
      

          $source=$gml_id{$source};
          $target=$gml_id{$target};
       


          print GML "edge [source ".$source." target ".$target." graphics [ fill \"\#".$line_color."\" width 4 ]]\n"; 

        }    
    }





    print GML "]\n";
    close GML;

}


sub initial_seed{
    my $class=shift;
    my $start_id;
foreach my $id(keys %{$class->{id_tn}}){
       $start_id=$id;
       last; 
     } 

    my ($I_id,$I_dist)=&max_dist($class,$start_id);
    my ($II_id,$II_dist)=&max_dist($class,$I_id);

   ## print $start_id."\t".$I_dist."\t".$I_id."\t".$II_dist."\t".$II_id."\n";
    return $class->{id_tn}->{$II_id};
}

sub pool_only{
    my ($class,$pool_ref)=@_;
    $class->{pool}=$pool_ref;
}


sub get_max_PD{
    my ($class,$limit,$print_ref)=@_;
    if (!($limit >= 1)){return;}
    my @array;
        
    my $i=1;

    while($i<=$limit){   
    my ($max_id,$max_dist,$remain_count)=$class->ident_max_PD();
   
    if(!$max_id){last;} 
    else{
         $class->grow_sub_tree($max_id);
         print $print_ref $class->{id_tn}->{$max_id}."\t".$max_dist."\n";
         $i++;
     
         if($remain_count <= 0) {last;}
       }
     }

    return;

   }

sub ident_max_PD{
    my $class=shift;
    my $max=-1e+100;
    my $max_id;
    my $remain_count=0;

    if($class->{pool}){      
    my @pool;
    foreach my $k(keys %{$class->{pool}}){
	if($class->{pool}->{$k}){
	    push (@pool,$k);
	}
    }
     
    foreach my $acc(@pool) {
	$remain_count++;
        my $id=$class->{tn_id}->{$acc};
    my $pd=&add_PD_single($class,$id);
    if($pd > $max){
	$max=$pd;
        $max_id=$id;
                  }
    }

    if($remain_count<=0){return;}
    $remain_count--;
    my $acc=$class->{id_tn}->{$max_id};   
    return ($max_id,$max,$remain_count);
}


	else{	  
    foreach my $id(keys %{$class->{remain_taxon}}) {
	$remain_count++;
    my $pd=&add_PD_single($class,$id);
    if($pd > $max){
	$max=$pd;
        $max_id=$id;
     }
    }
 
    if($remain_count<=0){return;}
    $remain_count--;
    return ($max_id,$max,$remain_count);
}

}




sub add_PD_single{
    my ($class,$id)=@_;
##if it is in sub tree return 0
    if($class->{sub_tree}->{$id}){
	return 0;
    }

##identify the ref path for this id
    my $path=$class->{ref_path}->{$id};
    my @path=split(/\s+/,$path);
##map sub tree nodes,see what is left
    my @add_path;

    while(@path) {
	my $step=pop @path;
        push(@add_path,$step);
	if($path=$class->{sub_tree}->{$step}){
	    last;
	}
    }

##calculate and return the left over PD
    my $path_len=scalar @add_path;
    $path_len-=2;
    my $PD;
    if($path_len>=0){
    for my $i(0..$path_len){
##	print $add_path[$i]."\t".$add_path[$i+1]."\t".$class->{matrix}->{$add_path[$i]}->{$add_path[$i+1]}."\n";
        $PD+=$class->{matrix}->{$add_path[$i]}->{$add_path[$i+1]}
    } 
    }
    else{
	die "Matrix BUG 001 \n";
    }
## return join(" ",@add_path);
    return $PD;    
}

sub grow_sub_tree{
   my ($class,$id)=@_;
       
	my @t=split(/\s+/,$class->{ref_path}->{$id});
            foreach my $t(@t){
	$class->{sub_tree}->{$t}=1;
        if($class->{remain_taxon}->{$t}){
	    delete $class->{remain_taxon}->{$t};}
             }
       	$class->{sub_tree}->{$id}=2;
   my $acc=$class->{id_tn}->{$id};
   if($class->{pool}){
   if($class->{pool}->{$acc}){
       $class->{pool}->{$acc}=0;
     }
   }
    return;
}


sub build_sub_tree{
    my ($class,@input_list)=@_;
    my $seed=$input_list[0];
    $class->ref_path_from_one_taxon($seed);

    foreach my $taxon(@input_list){
        my $id=$class->{tn_id}->{$taxon};
        if(!$id){
        die "$taxon from the list is not in the tree\n";
	}

	my @t=split(/\s+/,$class->{ref_path}->{$id});
            foreach my $t(@t){
	$class->{sub_tree}->{$t}=1;
        delete $class->{remain_taxon}->{$t};
        if($class->{pool}){
        my $acc=$class->{id_tn}->{$t};
        $class->{pool}->{$acc}=0; 
        }
             }
         }    
    return;
}


sub max_dist{
      my ($class,$taxon)=@_;
       if(!($taxon=~/T\d+/)) {
	   die "$taxon from the list is not in the tree\n";
       }

    my %include;
    my %current;
    my %dist;
    my $go_on=1;  

    $include{$taxon}=1;
    $current{$taxon}=1;    
    $dist{$taxon}=0;
    while($go_on){
        $go_on=0;

	foreach my $this_tn(keys %current){
	    if ($current{$this_tn}){
		foreach my $i(keys %{$class->{matrix}->{$this_tn}}){
		    if ($class->{matrix}->{$this_tn}->{$i} != -99999) {
			if ($include{$i}) {next;}
                        $dist{$i}=$dist{$this_tn}+$class->{matrix}->{$this_tn}->{$i};
                        $current{$i}=1;
                        $include{$i}=1;
                        if ($i=~/^N\d+/) {$go_on++;}  
                    }
		}
                }    
		$current{$this_tn}=0;             
        }
    }

      my $max_dist=-9999999;
      my $max_id;
      foreach my $id(keys %dist){
          if($id=~/^N\d+/){next;}
	  if($dist{$id}>=$max_dist){
	      $max_dist=$dist{$id};
              $max_id=$id;
	  }
      }

      undef %include;
      undef %current;
      undef %dist;      
          
      return ($max_id,$max_dist);

}

sub ref_path_from_one_taxon{
       my ($class,$taxon_str)=@_;
       my $taxon=$class->{tn_id}->{$taxon_str};
       if(!($taxon=~/T\d+/)) {
	   die "$taxon_str from the list is not in the tree\n";
       }


    my %include;
    my %current;
    my %path;

    my $go_on=1;  

    $include{$taxon}=1;
    $current{$taxon}=1;    
    $path{$taxon}=$taxon;

    while($go_on){

        $go_on=0;

	foreach my $this_tn(keys %current){
	    if ($current{$this_tn} == 1){
###problem here not all taxa are included
		foreach my $i(keys %{$class->{matrix}->{$this_tn}}){
		    if ($class->{matrix}->{$this_tn}->{$i} != -99999) {
			if ($include{$i}) {next;}
                        $path{$i}=$path{$this_tn}."\t"."$i";
                        $current{$i}=1;
                        $include{$i}=1;
                        if ($i=~/^N\d+/) {$go_on++;}  
                    }
		}
                }    
		
                delete $current{$this_tn};
                if($this_tn=~/^N\d+/){delete $path{$this_tn};}             
        }
    }

       
       foreach my $this_tn(keys %path){
        
	   if($this_tn=~/^N/){ 
           ##  print $path{$this_tn}."########\n";
             delete $path{$this_tn};
	   }
           else{$class->{remain_taxon}->{$this_tn}=1;}
       }
    
       $class->{ref_path}=\%path;
       return;
}

sub new{
    my ($class,$tree)=@_;   
    $tree=~s/\)\d+:/\):/g;
    $tree=~s/\|/__/g; ##take care of "|"
    $tree=~s/\"//g; ##take care of """
    $tree=~s/\'//g; ##take care of "'"
    $tree=~s/([A-Za-z][^\s,:;\(\)]+)\s+[A-Za-z][^\s,:;\(\)]+ *:/$1:/g; ##take care of space
  

    $tree=~s/\s+//g;
    $tree=~s/;$//;
   
    $tree=~s/\)\)/\):0.000\)/g;

    my $check=&node_number($tree);
    if (!$check) {return 0;}
    my $taxon_ref=&all_taxons($tree);

    ## setup taxon hash;
    my %tn_id;
    my %id_tn;
    my $i=0;
    my $id;
    my $tn;

    foreach $tn(@$taxon_ref){
        $i++;
	$id='T'.$i.'DONGYINGWU';

	$tn_id{$tn}=$id;
        $id_tn{$id}=$tn;

         $tn=~s/\|/\\\|/g;
         $tn=~s/\./\\\./g;
  
         $tree=~s/\b$tn\b/$id/;                
         $tn=$id;
	$tree=~s/$tn,/$tn:0.000,/;
	$tree=~s/$tn\)/$tn:0.000\)/;
    }

    $tree=~s/\),/\)0.000,/g;
     
         
  
    my $matrix;

 
    while ($tree=~/\(([NT]\d+DONGYINGWU:[^:,]+,[NT]\d+DONGYINGWU:[^:,\)]+(,[NT]\d+DONGYINGWU:[^:,\)]+)?)\)/) {

if($tree=~/\(([NT]\d+DONGYINGWU:[^:,]+,[NT]\d+DONGYINGWU:[^:,\)]+,[NT]\d+DONGYINGWU:[^:,\)]+)\)/){
my $this=$1;
my @t=split(/,/,$this);
my $that="(".$t[0].",".$t[1]."):0.000,".$t[2];
$tree=~s/$this/$that/;
next;
}


	$i++;
    my $this_node=$1;
	my $id='N'.$i.'DONGYINGWU';
       
        my $replace_node='\('.$this_node.'\)';

       

        $tree=~s/$replace_node/$id/;
       

	my @t=split(/,/,$this_node);

        foreach my $t(@t) {

	    my($child_id,$dis)=split(/:/,$t);
                  
            if(length($dis) >=1) { 
            $matrix->{$child_id}->{$id}=sprintf("%.10f",$dis);
            $matrix->{$id}->{$child_id}=sprintf("%.10f",$dis);
           }    
	}    
    }
  
    if(!($tree=~/^N\d+DONGYINGWU$/)){
        print "please make sure the taxa don't contain :;,()\"\'| \n";
        print $tree."\n";
	exit;
    }

############################

    my $self={};


  

       $self->{matrix}=$matrix;
       $self->{tn_id}=\%tn_id;
       $self->{id_tn}=\%id_tn;
   
    
    bless $self, $class;
    return $self;

}


### node_number($tree)
### return: 0 if the tree is not in the right format
###         non 0 for the number of nodes in the tree

sub node_number
{
    my $tree=shift;

    my $left=0;
    my $right=0;
   
    if (!( ($tree=~/^\(/) &&($tree=~/\)$/) )) {return 0;}
    
    while($tree=~/\(/g){
	$left++;
     
    }
    while($tree=~/\)/g){
	$right++;
  
    }

    if ($left == $right) {return $left;}
    else {return 0;}
}


### all_taxons($tree)
### return ref to an array of taxons
sub all_taxons{
    my $tree=shift;
    $tree=~s/,/ /g;
    $tree=~s/\(/ /g;
    $tree=~s/\)/ /g;

    $tree=~s/^\s+//;
    $tree=~s/\s+$//;

    my @taxons;

    my @temp=split(/\s+/,$tree);
  
  foreach my $taxon(@temp) {
    
	if ($taxon=~/^:/) {next;}
        else{($taxon)=split(/:/,$taxon);}

##        else {$taxon=~s/:(-*\d*\.*\d*$|\s*NaN$)//;
##        if($taxon=~/\w:/){
##         ($taxon)=split(/:/,$taxon); 
##       }
	if($taxon){
          push(@taxons,$taxon);
    }   
    }
    return \@taxons;
}


1;
