#!/usr/bin/env bash

for module in lbm_stencil lbm_propagate lbm_geometry lbm_domain global_coordinates data_view
do
    echo ""
    echo " ==== test module $module ===="
    ./test_module.sh $module | grep ': all_tests ='
    echo ""
done
