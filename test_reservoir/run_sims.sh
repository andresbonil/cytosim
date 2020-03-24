#!/usr/bin/env bash
set -e
current_directory=$(pwd)
for run_directory in "$@"
do
    if true
    then
        # Delete any previous data (if there was any)
        rm -f ${run_directory}/scan/run*/.DS_Store
        rm -f ${run_directory}/scan/run*/*.txt
        rm -f ${run_directory}/scan/run*/*.cmo
        rm -f ${run_directory}/scan/run*/*.cym
        rm -f ${run_directory}/scan/.DS_Store

        for j in ${run_directory}/scan/run*
        do
            if [ -d $j ]; then rmdir $j; fi
        done
        if [ -d ${run_directory}/scan ]; then rmdir ${run_directory}/scan; fi

        # Make the folder where all the simulation folders will be
        mkdir ${run_directory}/scan

        # Use preconfig.py to make a all the config files, you can change the number to the desired number of simulations

        python preconfig.py 20 ${run_directory}/config.cym.tpl ${run_directory}/scan

        # Use collect to put all the config files in subfolders with formatted names

        python collect.py ${run_directory}/scan/run%04i/config.cym ${run_directory}/scan/config????.cym

        python scan.py "${current_directory}/../bin/sim" nproc=6 ${run_directory}/scan/run????
    fi

    if true
    then

        python scan.py "${current_directory}/../bin/report single>single_report.txt" nproc=6 ${run_directory}/scan/run????
        python scan.py "python ${current_directory}/report.py" nproc=6 ${run_directory}/scan/run????

    fi

done