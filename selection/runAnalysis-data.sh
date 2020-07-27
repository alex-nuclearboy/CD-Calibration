#!/bin/bash
for X in `seq 45934 1 45999`; do
    input="${RUNS_DATA}/run_${X}"
    output="${OUTPUT_DATA}/DATA-SEC_calibration-run-${X}"
    if [ -e ${input} ]; then
        echo "run number $X..."
        if [ -e ${output}.root ];then
            echo "...was already done."
        else
            scriptname="analysis-run-$X.sh"
            if [ -e ${scriptname} ]; then
                echo "...already in process."
            else
                echo "... PROCESSING DATA $X ..."
                echo "#!/bin/bash" >> ${scriptname}
		echo "request_run ${X}" >> ${scriptname}
		echo >> ${scriptname}
                echo "cd $PWD" >> ${scriptname}
                echo "./main -mode raw -fin cluster:${input} -r ${X} -n ${output} -abort" >> ${scriptname}
                echo >> ${scriptname}
		#echo "release_run ${X}" >> ${scriptname}
                echo "rm -f $PWD/${scriptname}" >> ${scriptname}
                chmod u+x ${scriptname}
                qsub -q batch ${scriptname}
		sleep 2
                echo "...done."
            fi
        fi
    fi
done

