FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB metagenomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN gem install thinreports
ENV PATH /opt/conda/envs/metagenomics-profiling-1.0/bin:$PATH
RUN apt-get update && apt-get -y install build-essential
RUN mkdir /opt/samtools && \
	cd /opt/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar -xvf samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix=/opt/samtools/1.9 && \
	make install && \
	cd /opt/samtools && \
	rm -Rf samtools-1.9*
