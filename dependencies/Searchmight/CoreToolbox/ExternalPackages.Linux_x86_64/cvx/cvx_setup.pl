#!/usr/bin/perl

# just calls matlab to run the CVX setup

my $CMD;

open(CMD,"|matlab -nojvm > /dev/null");
print CMD "cvx_setup\n";
close CMD || warn "cvx_setup: matlab exited $?";
