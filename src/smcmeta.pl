#!/usr/bin/perl

# Sample command: ./sampler 100 -g ~/HOMD_16S_rRNA_RefSeq_V13.2.fasta -f 1000 -n 3 -s 0.2 0.9 -p -l 300 >profile
$cmd = "sampler " . join(" ", @sampler_pars);
# Simulator command: ./simulator -t 0.001 -r 0.001 -eG 1 1 -eg 2 1 0.1 -eM 3 0.5 -em 0.1 1 2 0.1 -eN 0.5 0.5 -en 0.4 0 0.01 -ej 0.2 0 1 < small_profile 2>debug.txt
# GTR model: ./simulator -t 0.001 -m GTR -R 1 0 0 0 0 1 < small_profile 2>debug.txt
$cmd = "| simulator " . join(" ", @simulator_pars);
`$cmd`;
