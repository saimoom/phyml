#!/usr/bin/perl
use File::Copy;
use Cwd;

$output_dir = $ARGV[0];
chdir($output_dir);
mkdir("./phyml/");
mkdir("./phyml/src");
mkdir("./phyml/doc");
mkdir("./phyml/bin");
mkdir("./phyml/examples");

# Put your path to the svn sources here
# $phyml_dir = "/home/lefort/homedir/home/ATGC/phyml/src/svn/phyml/";
$phyml_dir = "/Users/guindon/cvshome/phyml/";
chdir $phyml_dir;

open(SVNVERSION,"| svnversion > version");
close(SVNVERSION);


$phyml_files[0]  = "main";
$phyml_files[1]  = "utilities";
$phyml_files[2]  = "optimiz";
$phyml_files[3]  = "lk";
$phyml_files[4]  = "bionj";
$phyml_files[5]  = "models";
$phyml_files[6]  = "free";
$phyml_files[7]  = "options";
$phyml_files[8]  = "simu";
$phyml_files[9]  = "eigen";
$phyml_files[10] = "pars";
$phyml_files[11] = "alrt";
$phyml_files[12] = "interface";
$phyml_files[13] = "cl";
$phyml_files[14] = "spr";
$phyml_files[15] = "mc";
$phyml_files[16] = "m4";
$phyml_files[17] = "draw";
$phyml_files[18] = "rates";
$phyml_files[19] = "mg";
$phyml_files[20] = "mpi_boot";
$phyml_files[21] = "numeric";

chdir($phyml_dir); 
opendir(PHYML_DIR,$phyml_dir) || die;
@phyml_dir_files = readdir(PHYML_DIR);
closedir(PHYML_DIR);

foreach $filename (@phyml_dir_files)
  {
      foreach $source (@phyml_files)
      {
	  if($filename =~ m/$source.[c|h]$/)
	  {
	      copy $filename,$output_dir."/phyml/src/" || die "Could not copy file $filename.\n";
	  }
      }
  }

copy "Makefile.am",$output_dir."/phyml/src" || die "Could not copy file Makefile.am.\n";
copy "configure.ac",$output_dir."/phyml/src" || die "Could not copy file configure.ac.\n";
copy "INSTALL",$output_dir."/phyml/src" || die "Could not copy file INSTALL.\n";
copy "AUTHORS",$output_dir."/phyml/src" || die "Could not copy file AUTHORS.\n";
copy "README",$output_dir."/phyml/src" || die "Could not copy file README.\n";
copy "NEWS",$output_dir."/phyml/src" || die "Could not copy file NEWS.\n";
copy "COPYING",$output_dir."/phyml/src" || die "Could not copy file COPYING.\n";
copy "depcomp",$output_dir."/phyml/src" || die "Could not copy file depcomp.\n";
copy "install-sh",$output_dir."/phyml/src" || die "Could not copy file install-sh.\n";
copy "missing",$output_dir."/phyml/src" || die "Could not copy file missing.\n";
copy "ChangeLog",$output_dir."/phyml/src" || die "Could not copy file ChangeLog.\n";
copy "version",$output_dir."/phyml/src" || die "Could not copy file version.\n";

copy "phyml.bat",$output_dir."/phyml/bin" || die "Could not copy file phyml.bat.\n";
system("ps2pdf /Users/guindon/latex/phyml/phyml_manual.ps /Users/guindon/latex/phyml/phyml_manual.pdf");
copy("/Users/guindon/latex/phyml/phyml_manual.pdf",$output_dir."/phyml/doc") || die "Could not copy file phyml_manual.pdf.\n";
copy("/Users/guindon/data/phyml_example/nucleic",$output_dir."/phyml/examples/") || die "Could not copy example file 'nucleic'";
copy("/Users/guindon/data/phyml_example/proteic",$output_dir."/phyml/examples/") || die "Could not copy example file 'proteic'";


chdir($output_dir."/phyml/src");

# system("sed -e 's/#define VERSION.*/#define VERSION \"v3.0 '`head -n 1 version`'\"/' utilities.h > new_utilities.h");
# rename("new_utilities.h", "utilities.h");

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
# my @abbr = qw( January February March April May June July August September October November December );
my @abbr = qw( 01 02 03 04 05 06 07 08 09 10 11 12 );
$year = $year+1900;
$tarfile = "phyml_".$year.$abbr[$mon].$mday.".tar.gz";
chdir $output_dir;
system("tar -zcvf $tarfile ./phyml/ &> /dev/null");
print "<$mday>", "<$abbr[$mon]>", "<$year>","\n";
