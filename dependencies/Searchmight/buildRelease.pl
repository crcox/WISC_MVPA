#!/usr/bin/perl

#
# Prepares a release for source and all the architectures supported
#

use strict;

# defaults
my $architecture;
#my @architectures = qw/ Source /;
my @architectures = qw/ Darwin_i386 Linux_i686 Linux_x86_64 /;
#my @architectures = qw/ Linux_i686 /;
#my @architectures = qw/ Darwin_i386 /;
my $outputDir = "../SearchmightReleases";

my $numArgs = $#ARGV + 1;
my $version;
my @bits; my $bit; my $cmd;

if ($numArgs < 1) { die "syntax: prepareReleases.pl <version #>\n"; }

my $version = $ARGV[0];

#$architecture = @architectures[0];

my $cur = `pwd`;

foreach $architecture (@architectures) {

  # create output directory and subdirectories
  my $archDir = "SearchmightToolbox.$architecture.$version";
  my $outputDirVersion = "$outputDir/$archDir";
  mkdir($outputDirVersion);
  mkdir("$outputDirVersion/CoreToolbox");
  
  # copy all files (or the files they point to, if they are links) to output directory
  $cmd = "cp * $outputDirVersion";
  system($cmd);
  
  # copy the architecture specific directory
  $cmd = "cp -rL CoreToolbox/ExternalPackages.$architecture $outputDirVersion/CoreToolbox";
  system($cmd);
  
  # package into a tarball and remove directory
  
  chdir($outputDir);
  $cmd = "tar cvf $outputDir/searchmight.$architecture.$version.tar $archDir";
  system($cmd);
  $cmd = "gzip searchmight.$architecture.$version.tar";
  system($cmd);
  $cmd = "rm -rf $archDir";
 # system($cmd);
  $cmd = "ln -sf searchmight.$architecture.$version.tar.gz searchmight.$architecture.tar.gz";
  system($cmd);

  chdir('../SearchmightToolbox');
}

# update version

open OFD, "> $outputDir/VERSION";
print OFD "$version";
close OFD;
