FROM broadinstitute/gatk:4.1.4.1

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}

# Install miniconda environment
RUN conda update --all && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n mycodentifier python=3.6.8 pysam=0.15.2 nextflow=19.10.0 simplejson=3.17.0 pandas=1.0.1 matplotlib=3.1.2

# activate by bash mycodentifier env
RUN echo "source activate mycodentifier" > ~/.bashrc
ENV PATH /miniconda/envs/mycodentifier/bin:$PATH

SHELL ["/bin/bash", "-c"]

# add codebase to docker
ADD ./bin/snpit-master /workflow/bin/snpit-master

# adding snpit data from git and
RUN cd /workflow/bin/snpit-master && \
    python setup.py install && \
    python setup.py test

#run nextflow pipeline to verify if pipeline works
RUN cd /workflow && \
    mkdir -p /workflow/output
#    nextflow run /workflow/myco.nf --reads "/workflow/test/H37Rv_MB_R{1,2}.fastq.gz" --outDir "./output/TEST" --threads 8 --sampleName TEST
#
## remove work folder of nextflow test and conda environment
#RUN rm -rf /work && \
#    rm -rf /workflow/work && \
#    rm -rf /workflow/output/TEST && \
#    rm -rf /workflow/.git && \
#    rm -rf /workflow/conda && \
#    rm -rf /workflow/.nextflow && \
#    rm -rf /workflow/.nextflow.log && \
#    apt-get autoremove -y && \
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*