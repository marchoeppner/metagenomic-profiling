Bootstrap:docker
From:continuumio/anaconda

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the metagenomics-profile pipeline
    VERSION 1.0

%environment
    PATH=/opt/conda/envs/metagenomics-profiling-1.0/bin:$PATH
    export PATH

%files
    environment.yml /

%post

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
    gem install thinreports
    mkdir -p /ifs
    apt-get -y install procps

