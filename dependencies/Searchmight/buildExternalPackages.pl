#!/usr/bin/perl

use File::Basename;

#
# Creates an architecture specific version of all the external packages for this machine
#
# fpereira@princeton.edu
#

print "The toolbox requires some external packages that need to be compiled\n";
print "or installed for each specific combination of computer and operating system.\n";
print "I will proceed to build all I can automatically and give you instructions for\n";
print "what needs to be done manually at the end, if anything.\n";

my $location  = `pwd`;      chomp($location);
my $system    = `uname -s`; chomp($system);
my $processor = `uname -p`; chomp($processor);
my $release   = `uname -r`; chomp($release);

my $sourceDirectory = "ExternalPackages.template";
#my $targetDirectory = "ExternalPackages.$system\_$processor\_$release";
my $targetDirectory = "ExternalPackages.$system\_$processor";

system("rm -rf $targetDirectory");
system("cp -Rp $sourceDirectory $targetDirectory");
print "Created external package directory successfully!\n";

# create the MATLAB function that will set all the necessary paths
# for this combination of computer and OS

my $OFD;
open OFD, "> $targetDirectory/setupPaths.m";
print OFD "addpath $location/$targetDirectory/cvx/builtins\n";
print OFD "addpath $location/$targetDirectory/cvx/commands\n";
print OFD "addpath $location/$targetDirectory/cvx/functions\n";
print OFD "addpath $location/$targetDirectory/cvx/lib\n";
print OFD "addpath $location/$targetDirectory/cvx/structures\n";
print OFD "addpath $location/$targetDirectory/libsvm\n";
print OFD "addpath $location/$targetDirectory\n";
close OFD;
print "Created MATLAB path initialization successfully!\n";

# figure out where MATLAB is in this system
my $matlabFullPath; my $output;

$matlabFullPath = `readlink \`which matlab\``;
chomp($matlabFullPath);

my $discard; my $matlabPath;
($discard, $matlabPath, $discard) = fileparse($matlabFullPath);
chop($matlabPath);
($discard, $matlabPath, $discard) = fileparse($matlabPath);
chop($matlabPath);

print "Building packages:\n";

chdir("$targetDirectory");

# install CVX
chdir("cvx");
system("make");
print "- CVX - built successfully\n";
chdir("..");

# install libsvm (requires creating a new Makefile and running that)
my $OFD;
my $IFD;

chdir("libsvm");

open IFD, "<Makefile";
open OFD, ">Makefile.local";

while(<IFD>) {
    chomp;
    if (m/MATLABDIR \?=/) {
        print OFD "MATLABDIR ?= $matlabPath\n";
    } else {
        print OFD "$_\n";
    }
}

close IFD;
close OFD;

# and run it
system("cp Makefile.local Makefile");
system("make");

chdir("..");



# install anything else
chdir("..");

print "\nBuilt everything successfully!\n";
exit;

# instructions for packages that have to be configured manually

print "The rest has to be done manually. Please follow these steps:\n";
print "1) cd $targetDirectory\n";
print "2) cd libsvm\n";
print "3) edit Makefile so that MATLABDIR ?= <MATLAB main directory>\n";
print "   (use \"which matlab\" to find the matlab command, then\n";
print "   \"ls -l <matlab command path>\" to see its actual directory\n";
print "   and replace that in the Makefile\n";
print "4) make\n";
print "There are more details in ExternalPackages/README, if necessary\n";
