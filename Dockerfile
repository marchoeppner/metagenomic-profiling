FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB metagenomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/metagenomics-profiling-1.0/bin:$PATH
RUN apt-get update && \
	apt-get -y install build-essential ruby-full openssl libncurses5-dev zlib1g-dev xz-utils curl libcurl-dev bzip2 ca-certificates && \
	gem install thinreports
RUN mkdir /opt/samtools && \
	cd /opt/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar -xvf samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix=/opt/samtools/1.9 && \
	make install && \
	cd /opt/samtools && \
	rm -Rf samtools-1.9*

