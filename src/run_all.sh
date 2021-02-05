#!/usr/bin/env bash

for i in {1..33}; do
    echo
done

for path in \
    source \
    examples/multiphase_drop_impact \
    examples/multiphase_coalescence \
    examples/multiphase_coalescence_impact \
    examples/multiphase_drop_on_drop_impact \
    examples/multiphase_rising_bubble \
    examples/singlephase_couette_flow \
    examples/singlephase_poiseuille_flow_3D_channel \
    examples/singlephase_poiseuille_flow_plate \
    examples/singlephase_lid_driven_cavity \
; do

    echo
    echo
    echo
    echo $path
    echo

    ## remove old soft links
    rm -f ./openLBMflow_conf.c

    ## set new soft link
    ln -s ../$path/openLBMflow_conf.c .

    ## recompile
    make -j -B FEAFLAGS="-D TIME_TOTAL=200 -D TIME_SAVE=200"

    ## make sure we do not compare old files
    rm -rf output

    frame=0000200
    ## run and compare
    ./LBMflow $1 && \
    ../tools/numdiff  ./output/openLBMflow_$frame.vti \
            ../$path/reference/openLBMflow_$frame.vti   ## original source
#                    ./ref/$path/openLBMflow_$frame.vti  ## previous versions in src/ (see below)
    echo $path
    echo

    ## save a copy to compare for exact zero between code changes that should not affect anyting
#     mkdir -p ./ref/$path
#     mv ./output/openLBMflow_0000200.vti ./ref/$path

done


