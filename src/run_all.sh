#!/usr/bin/env bash

for i in {1..33}; do
    echo
done

#     examples/singlephase_couette_flow \
#     examples/singlephase_lid_driven_cavity \
for path in \
    source \
    examples/multiphase_drop_impact \
    examples/multiphase_coalescence \
    examples/multiphase_coalescence_impact \
    examples/multiphase_drop_on_drop_impact \
    examples/multiphase_rising_bubble \
    examples/singlephase_poiseuille_flow_3D_channel \
    examples/singlephase_poiseuille_flow_plate \
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
    
    ## run and compare
    ./LBMflow && \
    ../tools/numdiff  ./output/openLBMflow_0000200.vti \
            ../$path/reference/openLBMflow_0000200.vti   ## original source
#                    ./ref/$path/openLBMflow_0000200.vti  ## previous versions in src/
    echo $path
    echo

    ## save copy
#     mkdir -p ./ref/$path
#     mv ./output/openLBMflow_0000200.vti ./ref/$path

done


