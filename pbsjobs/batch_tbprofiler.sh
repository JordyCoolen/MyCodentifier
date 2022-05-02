#!/bin/bash

# first parameter: single file or full path of reads
if [[ -f "$1" ]]; then
    reads_path="$(readlink -f ${1})"
    reads_dir=${reads_path};
else
    reads_path="$(readlink -f ${1})"
    reads_dir=${reads_path}/*;
fi;

# second parameter: full path of output folder
out_path=${2:-"results/"}
out_dir="$(readlink -f ${out_path})"


for f in $reads_dir
do
    base=$(basename "$f" .fastq.gz)
        if [[ $f == *"_R"* ]]; then
            name="$(echo "$base" | awk '{ gsub(/_R+/, " " ); print $1; }')"
#            outdir=${2:-"./results/$name"}
            # a general glob for the two read pairs in nextflow
            reads="$(echo "$f" | awk '{ gsub(/_R+/, " " ); print $1; }')_R{1,2}.fastq.gz"
            # two separate read pairs for Tb-profiler
            read1="$(echo "$f" | awk '{ gsub(/_R+/, " " ); print $1; }')_R1.fastq.gz"
            read2="$(echo "$f" | awk '{ gsub(/_R+/, " " ); print $1; }')_R2.fastq.gz"

            if [[ ! " ${processed[@]} " =~ " ${reads} " ]]; then

                basepath="/home/gerwin/hseverin/project/tbprofiler_output/"

                basecmd="source activate tbprofiler\n\
                cd ${basepath}\n\
                mkdir -p ${out_dir}/${name}\n\
                { time tb-profiler profile -1 \"${read1}\" -2 \"${read2}\" -p \"${name}\" --dir \"${out_dir}/${name}\" -t 4 --txt 1> ${out_dir}/${name}/log.stderr ; } 2> ${out_dir}/${name}/time.txt"
                job_name="${name}$(date +%d_%m_%Y_%H_%M)"

                PBS="#!/bin/bash\n\

                #PBS -V\n\
                #PBS -N ${job_name}\n\
                #PBS -l nodes=1:ppn=12\n\
                #PBS -l walltime=2:00:00\n\
                #PBS -o ${job_name}.out\n\
                #PBS -e ${job_name}.err\n\
                cd \$PBS_O_WORKDIR\n\
                $basecmd"

                # Note that $PBS_O_WORKDIR is escaped ("\"). We don't want bash to evaluate this variable right now. Instead it will be evaluated when the command runs on the node.

                # Echo the string PBS to the function qsub, which submits it as a cluster job for you
                # A small delay is included to avoid overloading the submission process

                printf '\nQ-subbing job: %s \n' ${job_name}
                echo -e "$PBS" | awk '{$1=$1}1'

                echo -e ${PBS} | qsub
                sleep 0.5
                printf '\nsubmitted\n'

                processed+=( $reads );
        fi;
    else
            echo "reads name in $f does not have the form *_R{1,2}";
    fi;
done
