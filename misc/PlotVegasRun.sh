#! /bin/bash
Here=`pwd`

# comments:
# $1 file name with script, $2 max. number of iterations
# `(echo $1 | sed -e 's!/!\\\/!g')` replaces / with \// and . with \.
# check with echo `(echo $1 | sed -e 's!/!\\\/!g')`
# then sed replaces dummyFILE and dummyMAXIT


sed -e "s/dummyFILE/`(echo $1 | sed -e 's!/!\\\/!g')`/g" -e "s/dummyMAXIT/$2/g" $Here/misc/PlotVegasRun.tmpgpl | gnuplot
