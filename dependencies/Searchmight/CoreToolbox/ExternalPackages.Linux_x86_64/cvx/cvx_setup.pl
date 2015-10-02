#!/usr/bin/perl

# just calls matlab to run the CVX setup
my $MATLAB = $ARGV[0] . "/bin/matlab";
my $startmatlab = sprintf("|%s -nojvm > /dev/null", $MATLAB);
my $CMD;

open(CMD,$startmatlab);
print CMD "cvx_setup\n";
close CMD || warn "cvx_setup: matlab exited $?";
