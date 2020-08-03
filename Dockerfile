FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y build-essential wget git autoconf
RUN apt-get install -y libgomp1 ncbi-blast+ hmmer cd-hit
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels r

RUN conda create -n LTR_retriever
RUN conda install -y -c conda-forge perl perl-text-soundex
RUN conda install -y -c bioconda cd-hit repeatmasker
## Install Repeat Masker libraries steps included courtesy @wangshun1211 for contributions
RUN cd /opt && wget https://github.com/lfaino/LoReAn/raw/noIPRS/third_party/software/RepeatMasker.Libraries.tar.gz && \
    tar -zxf RepeatMasker.Libraries.tar.gz && rm -rf /opt/conda/share/RepeatMasker/Libraries && \
    mv ./Libraries /opt/conda/share/RepeatMasker/. && chmod -R 755 /opt/conda/share/RepeatMasker/Libraries && \
    rm -f RepeatMasker.Libraries.tar.gz 

RUN git clone https://github.com/oushujun/LTR_retriever.git
RUN echo "source activate LTR_retriever" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
ENV PATH /LTR_retriever/bin:$PATH
ENV PATH /LTR_retriever:$PATH

ENTRYPOINT [ "/LTR_retriever/LTR_retriever" ]
