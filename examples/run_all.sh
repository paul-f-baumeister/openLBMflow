#!/bin/bash

src=../../source/openLBMflow.c

for d in *                                                            ## all cases
# for d in singlephase_*                                                ## only singlephase_ cases
# for d in multiphase_*                                                 ## only multiphase_ cases
do
  echo $d
  if [ -d $d ]; then ## must be a directory
    cd $d
      echo 
      echo "====================================================================="
      echo $d
      echo "====================================================================="
      pwd                                                               ## show the current directory
      echo

      rm -f openLBMflow.exe                                             ## clean up from previous runs
      rm -f ./openLBMflow.c                                             ## remove any old soft link
      ln -s $src ./openLBMflow.c                                        ## create a soft link to the source
      ## compile and execute
      gcc -O2 ./openLBMflow.c -lm -o openLBMflow.exe
      time  ./openLBMflow.exe                           

      rm -f openLBMflow.exe                                             ## clear executable
      rm ./openLBMflow.c                                                ## remove soft linke to the source
      mv output reference                                               ## rename output directory

    cd ..
  fi
done
