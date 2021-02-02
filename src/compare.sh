#!/usr/bin/env bash

reference=../source/output

../tools/numdiff  ./output/openLBMflow_0000200.vti \
                $reference/openLBMflow_0000200.vti
