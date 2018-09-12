#!/bin/bash

dir=$1
sed "s/dir/$dir/" runstats.sl > tmp.sl
DATE=$dir OUTPUT=$dir.stats.fits sbatch tmp.sl
\rm tmp.sl
