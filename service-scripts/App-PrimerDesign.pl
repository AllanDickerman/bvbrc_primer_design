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
use primer3plusFunctions;
use settings;
use HtmlFunctions;

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

    #my @command = ("primer3_core", "--format_output", "$p3params_file");
    #print "run command: ", join("\t", @command), "\n";
    my $cwd = getcwd();
    chdir("/homes/allan/git/dev_container/modules/bvbrc_primer_design/primer3plus-2.5.0/cgi-bin"); 
    my %resultsHash = {};
    $resultsHash{"TEST1_KEY"}="VALUE";
    my %defaultHash = {};
    my %defaultSettings = getDefaultSettings();
    #for my $key (keys %defaultSettings) {
    #    print STDERR "DefSet: $key $defaultSettings{$key}\n";
    #}
    #print STDERR "\n";
    # Add missing parameters from the defaultSettings
    my %completeParameters = constructCombinedHash( %defaultSettings, %$params );

    # check for parameters which make no sense and correct them
    checkParameters(%completeParameters);

    runPrimer3(\%completeParameters, \%defaultSettings, \%resultsHash);
    #my ($out, $err) = run_cmd(\@command);
    #print STDERR "STDOUT:\n$out\n";
    #print STDERR "STDERR:\n$err\n";
    print STDERR "Now showing results from runPrimer3\n";
    for my $key (keys %resultsHash) {
        print "results: $key\tvalue=$resultsHash{$key}\n";
    }
    print STDERR "Now show messages:\n";
    my @messages = getMessages();
    for my $msg (@messages) {
        print STDERR "Msg: $msg\n"
    }
	print "Content-type: text/html\n\n";
	my $html = mainResultsHTML( \%completeParameters, \%resultsHash ), "\n";

    my $html_file = "$tmpdir/$params->{output_file}_report.html";
    open F, ">$html_file";
    print F $html;
    close F;
    writeStatistics("primer3plus_run_primer3");
    chdir($cwd);
    exit();

    write_report($html_file);
    push @outputs, [$html_file, "html"];

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
    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");
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
