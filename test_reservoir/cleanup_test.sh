#!/usr/bin/env bash

for run_directory in "$@"
do
    # Delete any previous data (if there was any)
    rm -f ${run_directory}/scan/run*/.DS_Store
    rm -f ${run_directory}/scan/run*/*.cym
    rm -f ${run_directory}/scan/run*/*.cmo
    rm -f ${run_directory}/scan/run*/*.txt
    rm -f ${run_directory}/scan/.DS_Store

    for j in ${run_directory}/scan/run*
        do
            if [ -d $j ]; then rmdir $j; fi
        done
        if [ -d ${run_directory}/scan ]; then rmdir ${run_directory}/scan; fi

done