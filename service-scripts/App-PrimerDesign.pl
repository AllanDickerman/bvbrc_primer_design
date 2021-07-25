#
# The PrimerDesign application.
use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use JSON;
use File::Copy ('copy', 'move');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use URI::Escape;
#use primer3plusFunctions;
#use settings;
#use HtmlFunctions;

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

our $debug = 0;
$debug = $ENV{"PrimerDesignDebug"} if exists $ENV{"PrimerDesignDebug"};
print STDERR "debug = $debug\n" if $debug;
print STDERR "args = ", join("\n", @ARGV), "\n" if $debug;
our @analysis_step => ();# collect info on sequence of analysis steps
our @step_stack => (); # for nesting of child steps within parent steps

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#my $data_url = "http://www.alpha.patricbrc.org/api";
print STDERR "data_url=\n$data_url\n" if $debug;

my $script = Bio::KBase::AppService::AppScript->new(\&design_primers, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";
    my $pf = {
	cpu => 4,
	memory => "32G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub design_primers {
    my ($app, $app_def, $raw_params, $params) = @_;

    my $time1 = `date`;
    print "Proc DesignPrimers ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my @outputs; # array of tuples of (filename, filetype)
    my $api = P3DataAPI->new;
    my $tmpdir = File::Temp->newdir( "/tmp/PrimerDesign_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
   
    run("echo $tmpdir && ls -ltr $tmpdir");

    my $p3params_file = "$params->{output_file}_P3_parameters.txt";
    open F, ">$tmpdir/$p3params_file";
    print STDERR "params->{parameters} = $params\n";
    for my $param (keys %{$params}) {
        if ($param !~ /^output/) {
            my $value = $params->{$param};
            print F "${param}=$value\n";
        }
    }
    print F "=\n";
    close F;

    my $cwd = getcwd();
    #chdir("/homes/allan/git/dev_container/modules/bvbrc_primer_design/primer3plus-2.5.0/cgi-bin"); 
    chdir($tmpdir);
    #my %defaultSettings = getDefaultSettings();
    #for my $key (keys %defaultSettings) {
    #    print STDERR "DefSet: $key $defaultSettings{$key}\n";
    #}
    #print STDERR "\n";
    # Add missing parameters from the defaultSettings
    #my %completeParameters = constructCombinedHash( %defaultSettings, %$params );

    # check for parameters which make no sense and correct them
    #checkParameters(%completeParameters);

    #runPrimer3(\%completeParameters, \%defaultSettings, \%resultsHash);
    my $primer3_output_file = "$params->{output_file}_raw_output.txt";
    #my @command = ("primer3_core --output $primer3_output_file $p3params_file");
    my $command = "primer3_core --output $primer3_output_file";
    print "run command: $command\n";
    open PROC, "|$command";
    for my $param (keys %{$params}) {
        if ($param !~ /^output/) {
            my $value = $params->{$param};
            print PROC "${param}=$value\n";
            print "PROC ${param}=$value\n";
        }
    }
    print PROC "=\n";
    close PROC;
    #my ($out, $err) = run_cmd(\@command);
    #print STDERR "STDOUT:\n$out\n";
    #print STDERR "STDERR:\n$err\n";
    my %resultsHash = {};
    #$resultsHash{"TEST1_KEY"}="VALUE";
    print STDERR "Now showing results from runPrimer3\n";
    open F, $primer3_output_file;
    while (<F>)
    {
        chomp;
        my ($key, $val) = split("=", $_, '2');
        $resultsHash{$key} = $val;
        print "results: $key\tvalue=$resultsHash{$key}\n";
    }
    print STDERR "Now show messages:\n";
    #my @messages = getMessages();
    #for my $msg (@messages) {
    #    print "msg: $msg\n";
    #}
    my $pair_count = 0;
    # Figure out if any primers were returned and
    # write a helping page if no primers are returned  
    if (defined $resultsHash{"PRIMER_PAIR_NUM_RETURNED"}){
        $pair_count = $resultsHash{"PRIMER_PAIR_NUM_RETURNED"};
    }
    my $html = "<html>\n";
    
    if ($pair_count > 0) {
        my %html_param;
        $html_param{right_primer_color} = "rgb(250,240,75)";
        $html_param{left_primer_color} = "#ccccff";
        $html_param{row_alt_color} = "#dddddd";
        $html_param{row_highlight_color} = "#9999ff";
        $html .= generate_javascript(\%html_param);

        $html .= generate_html_all_pairs_view(\%resultsHash, \%html_param);	
    }
    else {
        $html = "Problem: no primer pairs found.\n";
        if ($resultsHash{ERROR}) {
            $html .= "<p>Error: $resultsHash{ERROR}\n";
        }
    }
    $html .= "</html>\n";
    my $html_file = "$tmpdir/$params->{output_file}_report.html";
    open F, ">$html_file";
    print F $html;
    #print F writeStatistics("primer3plus_run_primer3");
    close F;

    push @outputs, [$html_file, "html"];
    push @outputs, [$primer3_output_file, "txt"];

    print STDERR '\@outputs = '. Dumper(\@outputs);
    my $output_folder = $app->result_folder();
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n";
        if (0) {
           $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
               (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
               $global_token);
        }
        else {
            my $ext = $1 if $ofile =~ /.*\.(\S+)$/;
            my @cmd = ("p3-cp", "-f", "-m", "${ext}=$type", $ofile, "ws:" . $app->result_folder);
            print STDERR "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        }
    }
    chdir($cwd);
    my $time2 = `date`;
    #print STDERR ("Start: $time1"."End:   $time2"r. "$tmpdir/DONE\n");
}

sub generate_javascript {
    my $html_param = shift;
    return qq(
<script type="text/javascript">

let cur_pair_index = 0;
let right_primer_color="$html_param->{right_primer_color}";
let left_primer_color="$html_param->{left_primer_color}";
let alternate_row_colors = ["transparent", "$html_param->{row_alt_color}"];
let row_highlight_color="$html_param->{row_highlight_color}";

function update_current_pair() {
new_index = document.getElementById("select_pair_control").value;

old_row_color = alternate_row_colors[cur_pair_index % 2]
old_row = document.getElementById("PRIMER_LEFT_"+cur_pair_index);
old_row.style.backgroundColor=old_row_color;
old_row.children[4].style.backgroundColor="transparent";

old_row = document.getElementById("PRIMER_RIGHT_"+cur_pair_index);	
old_row.style.backgroundColor=old_row_color;
old_row.children[3].style.backgroundColor="transparent";

new_row = document.getElementById("PRIMER_LEFT_"+new_index);
new_row.style.backgroundColor=row_highlight_color;
new_row.children[4].style.backgroundColor=left_primer_color;

new_row = document.getElementById("PRIMER_RIGHT_"+new_index);
new_row.style.backgroundColor=row_highlight_color;
new_row.children[3].style.backgroundColor=right_primer_color;

document.getElementById("sequence_for_pair_"+cur_pair_index).style.display="none";
document.getElementById("sequence_for_pair_"+new_index).style.display="block";

document.getElementById("primer_pair_id").innerHTML=new_index;
document.getElementById("pair_data_"+cur_pair_index).style.backgroundColor="transparent";
document.getElementById("pair_data_"+new_index).style.backgroundColor=row_highlight_color;

cur_pair_index = new_index;
}

</script>
)
}

sub generate_html_all_pairs_view {
    my $results = shift;
    my $html_param = shift;
    my $pair_count = $results->{PRIMER_PAIR_NUM_RETURNED};
    my $html .= "<div id=\"primer3_results_container\" onload='update_current_pair(0)'>\n";
    $html .= "<b>Show Primer Pair:</b> <select id=\"select_pair_control\" onchange='update_current_pair(this.value)'>\n";
    for (my $i=0; $i < $pair_count; $i++) {
        $html .= "<option value=\"$i\">$i</option>\n";
    }
    $html .= "</select><p>\n";

    $html .= "<table class='bvbrc_primer_table' id='PRIMER_PAIR_TABLE'>\n";
    $html .= "<tr><th>Sel</th><th>Dir</th><th>Start</th><th>Length</th><th>Sequence</th><th>Tm</th><th>GC</th><th>Any</th><th>End</th><th>3' Stab</th><th>Penalty</th></tr>\n";
    for my $primer_index (0 .. $pair_count-1) {
        my $row_color = $primer_index % 2 ? $html_param->{row_alt_color} : "transparent";
        $row_color = $html_param->{row_highlight_color} unless $primer_index;
        for my $dir ("LEFT", "RIGHT") {
            my $key_prefix = "PRIMER_${dir}_$primer_index";
            $html .= "<tr id=\"$key_prefix\" style=\"background-color:$row_color\">";
            if ($dir eq 'LEFT') {
                $html .= "<td rowspan='2'>$primer_index</td>"; 
            }
            $html .= "<td>$dir</td>";
            my ($start, $length) = split(",", $results->{$key_prefix});
            $html .= "<td>$start</td>";
            $html .= "<td>$length</td>";
            my $seq_bg_color = "transparent";
            if ($primer_index == 0) {
                $seq_bg_color = $dir eq "LEFT" ? $html_param->{left_primer_color} : $html_param->{right_primer_color};
            }
            $html .= "<td style='background-color:$seq_bg_color'>" . $results->{$key_prefix . "_SEQUENCE"} . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_TM"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_GC_PERCENT"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_ANY"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_END"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_END_STABILITY"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_PENALTY"}) . "</td>";
            $html .= "</tr>\n";
        }
    }
    $html .= "</table>\n";

    $html .= "<p>Statistics for primer pair <span id='primer_pair_id'>0</span>\n";
    $html .= "<table  onload='update_current_pair()'>\n";
    $html .= "<tr><th>product Size</th><th>End Coml.</th><th>Any Compl.</th><th>Penalty</th></tr>\n";
    for my $pair_index (0 .. $pair_count-1) {
        my $row_color = $pair_index % 2 ? $html_param->{row_alt_color} : "transparent";
        $html .= "<tr id='pair_data_$pair_index' style='background-color:$row_color'>";
        my $prefix = "PRIMER_PAIR_$pair_index";
        $html .= "<td>".$results->{$prefix . "_PRODUCT_SIZE"}."</td>";
        $html .= "<td>".$results->{$prefix . "_COMPL_END_TH"}."</td>";
        $html .= "<td>".$results->{$prefix . "_COMPL_ANY_TH"}."</td>";
        $html .= "<td>".$results->{$prefix . "_PENALTY"}."</td></tr>\n";
    }
    $html .= "</table>\n";
    $html .= "<p>Primers in template sequence:<div id=\"sequence_container\" style=\"font-family:Courier,monospace\">\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
        $html .= create_primer_pair_on_sequence_html($results, $index, $index == 0); 
    } 
    $html .= "</div></div>\n"; 
    return $html;
}


sub generate_html_one_pair_view {
    my $results = shift;
    my $pair_count = $results->{PRIMER_PAIR_NUM_RETURNED};
    my $html .= "<b>Show Primer Pair:</b> <select id=\"select_pair_control\" onchange=\"update_current_pair(this.value)\">\n";
    for (my $i=0; $i < $pair_count; $i++) {
        $html .= "<option value=\"$i\">$i</option>\n";
    }
    $html .= "</select><p>\n";
    #$html .= "<span id=\"primer_pair_label\" style=\"font-weight: bold\">Showing Primer Pair 0</span>\n";
    $html .= "<div id=\"primer_pair_container\">\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
      $html .= create_primer_pair_table_html($results, $index, $index == 0); 
    } 
    $html .= "</div>\n"; 

    $html .= "<p>Primers in template sequence:<div id=\"sequence_container\" style=\"font-family:Courier,monospace\">\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
      $html .= create_primer_pair_on_sequence_html($results, $index, $index == 0); 
    } 
    $html .= "</div>\n"; 
    return $html
}

    

sub create_primer_pair_table_html {
  my ($results, $primer_index, $visibility) = @_; 

  my $tableHTML .= "<table class='bvbrc_primer_table' id='PRIMER_PAIR_TABLE_$primer_index' style=\"display:";
  $tableHTML .= $primer_index eq 0 ? "block" : "none";  # first one is visible to start (changed by javascript on user interaction)
  $tableHTML .= "\">\n";
  $tableHTML .= "<tr><th>Dir</th><th>Start</th><th>Length</th><th>Sequence</th><th>Tm</th><th>GC</th><th>Any</th><th>End</th><th>3' Stab</th><th>Penalty</th></tr>\n";
  for my $dir ("LEFT", "RIGHT") {
      my $key_prefix = "PRIMER_${dir}_$primer_index";
      $tableHTML .= "<tr>";
      $tableHTML .= "<td>$dir</td>";
      my ($start, $length) = split(",", $results->{$key_prefix});
      $tableHTML .= "<td>$start</td>";
      $tableHTML .= "<td>$length</td>";
      my $color = $dir eq "LEFT" ? "#ccccff" : "rgb(250,240,75)";
      $tableHTML .= "<td style=\"background-color:$color\">" . $results->{$key_prefix . "_SEQUENCE"} . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_TM"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_GC_PERCENT"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_ANY"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_END"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_END_STABILITY"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_PENALTY"}) . "</td>";
      $tableHTML .= "</tr>\n";
  }
  my $product_length = $results->{"PRIMER_PAIR_${primer_index}_PRODUCT_SIZE"};
  $tableHTML .= "<tr><td colspan=\"10\"><span style=\"font-weight:bold\">Product Length:</span> $product_length</td></tr>\n";
  $tableHTML .= "</table>\n";
  return $tableHTML;
}

###################################################
# divHTMLsequence: Prints out the sequence nicely # 
###################################################
sub create_primer_pair_on_sequence_html {
  my ($results, $sequence, $seqLength, $firstBase);
  my ($base, $count, $preCount, $postCount, $printBase) ;
  my ($format, $baseFormat, $preFormat, $pair_to_show);
  my (@targets, $region, $run, $counter, $madeRegion);
  my $tableHTML;
   
  $results = shift;
  $pair_to_show = shift;
  my $is_visible = shift;

  $sequence = $results->{"SEQUENCE_TEMPLATE"};
  $format = $sequence;
  $format =~ s/\w/N/g;

  $seqLength = length ($sequence);
  $firstBase = $results->{"PRIMER_FIRST_BASE_INDEX"};
  
  if ((defined $results->{"SEQUENCE_TEMPLATE"})and ($results->{"SEQUENCE_TEMPLATE"} ne "")) {

  if (defined ($results->{"SEQUENCE_EXCLUDED_REGION"}) and (($results->{"SEQUENCE_EXCLUDED_REGION"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_EXCLUDED_REGION"};
      foreach $region (@targets) {
          $format = addRegion($format,$region,$firstBase,"E");
      }
  }
  if (defined ($results->{"SEQUENCE_TARGET"}) and (($results->{"SEQUENCE_TARGET"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_TARGET"};
      foreach $region (@targets) {
          $format = addRegion($format,$region,$firstBase,"T");
      }
  }
  if (defined ($results->{"SEQUENCE_INCLUDED_REGION"}) and (($results->{"SEQUENCE_INCLUDED_REGION"}) ne "")) {
      $format = addRegion($format,$results->{"SEQUENCE_INCLUDED_REGION"},$firstBase,"I");
  }
  
  # Add only the specified primer pair e.g. Detection
  if ($pair_to_show >= 0) {
      if (defined ($results->{"PRIMER_LEFT_" . $pair_to_show}) and (($results->{"PRIMER_LEFT_" . $pair_to_show}) ne "")) {
           $format = addRegion($format,$results->{"PRIMER_LEFT_" . $pair_to_show},$firstBase,"F");
      }
      if (defined ($results->{"PRIMER_INTERNAL_" . $pair_to_show}) and (($results->{"PRIMER_INTERNAL_" . $pair_to_show}) ne "")) {
     	   $format = addRegion($format,$results->{"PRIMER_INTERNAL_" . $pair_to_show},$firstBase,"O");
      }
      if (defined ($results->{"PRIMER_RIGHT_" . $pair_to_show}) and (($results->{"PRIMER_RIGHT_" . $pair_to_show}) ne "")) {
           $format = addRegion($format,$results->{"PRIMER_RIGHT_" . $pair_to_show},$firstBase,"R");
      }
  }
  # Mark no primers on the sequence e.g. Primer list
  elsif ($pair_to_show eq -1) {
      
  }
  if (defined ($results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"}) and (($results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"};
      foreach $region (@targets) {
          $madeRegion = $region - 1;
          $madeRegion .= ",2";
          $format = addRegion($format,$madeRegion,$firstBase,"-");
      }
  }

  ## Handy for testing:
  #  $sequence = $format;

  my $display_style = $is_visible ? "block" : "none";
  $tableHTML = "<table class=\"primer3plus_table_no_border\" id=\"sequence_for_pair_$pair_to_show\" style=\"display: $display_style\">";
  $tableHTML .= qq{
     <colgroup>
       <col width="13%" style="text-align: right;">
       <col width="87%">
     </colgroup>
};
  $preFormat = "N";
  for (my $i=0; $i<$seqLength; $i++) {
     $count = $i;
     $preCount = $i - 1;
     $postCount = $i + 1;
     $base = substr($sequence,$i,1);
     $baseFormat = substr($format,$i,1);

     if (($count % 50) eq 0) {
         $printBase = $i + $firstBase;
         $tableHTML .= qq{     <tr>
       <td>$printBase&nbsp;&nbsp;</td>
       <td>};
     }

     if ($preFormat ne $baseFormat) {
         if ($preFormat ne "J") {
             $tableHTML .= qq{</a>};
         }
         if ($baseFormat eq "N") {
             $tableHTML .= qq{<a>};
         }
         if ($baseFormat eq "E") {
             $tableHTML .= qq{<a class="primer3plus_excluded_region">};
         }
         if ($baseFormat eq "T") {
             $tableHTML .= qq{<a class="primer3plus_target">};
         }
         if ($baseFormat eq "I") {
             $tableHTML .= qq{<a class="primer3plus_included_region">};
         }
         if ($baseFormat eq "F") {
             $tableHTML .= qq{<a style="background-color:#ccccff; color:rgb(0,0,0);">};
         }
         if ($baseFormat eq "O") {
             $tableHTML .= qq{<a class="primer3plus_internal_oligo">};
         }
         if ($baseFormat eq "R") {
             $tableHTML .= qq{<a style="background-color:rgb(250,240,75); color:rgb(0,0,0);">};
         }
         if ($baseFormat eq "B") {
             $tableHTML .= qq{<a class="primer3plus_left_right_primer">};
         }
         if ($baseFormat eq "-") {
             $tableHTML .= qq{<a class="primer3plus_primer_overlap_pos">};
         }

     }

     if ((($count % 10) eq 0) and !(($count % 50) eq 0)) {
         $tableHTML .= qq{&nbsp;&nbsp;};
     }

     $tableHTML .= qq{$base};

     if (($postCount % 50) eq 0) {
         $tableHTML .= qq{</a></td>
     </tr>
};
         $baseFormat = "J";
     }
     $preFormat = $baseFormat;
  }

  if (($postCount % 50) ne 0) {
      $tableHTML .= qq{</a></td>
     </tr>
};
  }

  $tableHTML .= "  </table>\n";
}

  return $tableHTML;
}

sub addRegion {
  my ($formatString, $region, $firstBase, $letter);
  my ($regionStart, $regionLength, $regionEnd);
  my ($stingStart, $stringRegion, $stringEnd, $stringLength);
  $formatString = shift;
  $region = shift;
  $firstBase = shift;
  $letter = shift;
  ($regionStart, $regionLength) = split "," , $region;
  $regionStart =~ s/\s//g;
  $regionLength =~ s/\s//g;
  $regionStart = $regionStart - $firstBase;
  if ($regionStart < 0) {
      $regionStart = 0;
  }
  if ($letter eq "R") {
  $regionStart = $regionStart - $regionLength + 1;  
  }
  $regionEnd = $regionStart + $regionLength;
  $stringLength = length $formatString;
  if ($regionEnd > $stringLength) {
      $regionEnd = $stringLength;
  }
  $stingStart = substr($formatString,0,$regionStart);
  $stringRegion = substr($formatString,$regionStart,$regionLength);
  $stringEnd = substr($formatString,$regionEnd,$stringLength);
  if (($letter ne "F") and ($letter ne "R")) {
      $stringRegion =~ s/\w/$letter/g;
  }
  if ($letter eq "F") {
      $stringRegion =~ tr/NETIFROB/FFFFFBFB/;
  }
  if  ($letter eq "R") {
      $stringRegion =~ tr/NETIFROB/RRRRBRRB/;
  }
  $formatString = $stingStart;
  $formatString .= $stringRegion;
  $formatString .= $stringEnd;

  return $formatString;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}
