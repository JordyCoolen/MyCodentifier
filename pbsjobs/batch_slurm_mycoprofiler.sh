#!/bin/bash

# USAGE
# bash batch_nextflow.sh <inputpath> <outputpath>

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

for f in ${reads_dir}
do
    base=$(basename "$f" .fastq.gz)
        if [[ $f == *"_R"* ]]; then
            name="$(echo "$base" | awk '{ gsub(/_R+/, " " ); print $1; }')"
            # a general glob for the two read pairs in nextflow
            reads="$(echo "$name" | awk '{ gsub(/_R+/, " " ); print $1; }')_R{1,2}.fastq.gz"
            # two separate read pairs for Tb-profiler
            read1="$(echo "$f" | awk '{ gsub(/_R+/, " " ); print $1; }')_R1.fastq.gz"
            read2="$(echo "$f" | awk '{ gsub(/_R+/, " " ); print $1; }')_R2.fastq.gz"

            if [[ ! " ${processed[@]} " =~ " ${reads} " ]]; then
                basecmd="mkdir -p ${out_dir}/${name}\n\
                (git log -1 && echo -e "\n") > ${out_dir}/${name}/log.stderr\n\
                { time singularity exec -B ${reads_path}:/data,/ifs/software/external/mmb/mycoprofiler:/workflow,/ifs/software/external/mmb/tools/centrifuge/db/:/workflow/db/centrifuge_WGS,/ifs/software/external/mmb/tools/centrifuge/db/:/workflow/db/centrifuge_cont /ifs/software/external/mmb/singularity_images/mycoprofiler_v1_72.sif nextflow run /workflow/myco.nf --threads 16 --reads /data/\"$reads\" --sampleName \"$name\" --outDir \"$out_dir\" --subsampling true --snpeff_base snpEff/ 1 >> ${out_dir}/${name}/log.stderr ; } 2> ${out_dir}/${name}/time.txt &&\n\
                out_dir=${out_dir} &&\n\
                name=${name} &&\n\
                cat \${out_dir}/\${name}/*/.command.log > \${out_dir}/\${name}/\${name}.cmd.log &&\n\
                idline=\$(grep "Launching.*[.*]" \${out_dir}/\${name}/log.stderr) &&\n\
                logfile=\$(grep -HF \"\${idline}\" .nextflow.log* | cut -d : -f 1) &&\n\
                mv \$logfile \${out_dir}/\${name}/\${name}.nf.log"

                job_name="${name}$(date +%d_%m_%Y_%H_%M)"

                PBS="#!/bin/bash\n\
                #SBATCH -J ${job_name}\n\
                #SBATCH --ntasks=1\n\
                #SBATCH --cpus-per-task=16\n\
                #SBATCH --mem-per-cpu=16100M\n\
                #SBATCH —output=${job_name}.out\n\
                #SBATCH —error=${job_name}.err\n\
                #SBATCH --partition=research\n\
                cd /ifs/software/external/mmb/singularity_images\n\
                module load gcc\n\
                module load bioinf/singularity/3.5.3\n\
                $basecmd"

                # Note that $PBS_O_WORKDIR is escaped ("\"). We don't want bash to evaluate this variable right now. Instead it will be evaluated when the command runs on the node.

                # Echo the string PBS to the function qsub, which submits it as a cluster job
                # A small delay is included to avoid overloading the submission process

                #echo -e "${PBS}"
                printf '\nQ-subbing job: %s \n' ${job_name}
                echo -e "$PBS" | awk '{$1=$1}1'

                echo -e ${PBS} | sbatch
                sleep 0.5
                printf '\nsubmitted\n'

                processed+=( $reads );
        fi;
    else
            echo "reads name in $f does not have the form *_R{1,2}";
    fi;
done