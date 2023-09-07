# 11.3.1-base-ubuntu20.04
FROM nvidia/cuda:11.8.0-base-ubuntu22.04 as base 
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG C.UTF-8


RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get -qq update \
    && apt-get -qq install -y gnupg ca-certificates wget git \
    && apt-get -qq clean \
    && rm -rf /var/lib/apt/lists/*


RUN  export DEBIAN_FRONTEND=noninteractive \
    && wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb \
    && dpkg -i cuda-keyring_1.1-1_all.deb \
    && apt-get -qq update \ 
    && apt-get install -y --no-install-recommends \
        cuda-command-line-tools-11-8 \
        cuda-cudart-dev-11-8 \ 
        cuda-nvcc-11-8 \
        cuda-cupti-11-8 \
        cuda-nvprune-11-8 \
        cuda-libraries-11-8 \
        cuda-nvrtc-11-8 \ 
        libcufft-11-8 \
        libcurand-11-8 \
        libcusolver-dev-11-8 \ 
        libcusparse-dev-11-8 \
        libcublas-dev-11-8 \
        libcudnn8=8.6.0.163-1+cuda11.8 \
        libnvinfer-plugin8=8.6.1.6-1+cuda11.8 \
        libnvinfer8=8.6.1.6-1+cuda11.8 \
        build-essential \
        pkg-config \
        curl \
        software-properties-common \
        unzip \
    && apt-get -qq clean \
    && rm -rf /var/lib/apt/lists/* 

RUN find /usr/local/cuda-*/lib*/ -type f -name 'lib*_static.a' -not -name 'libcudart_static.a' -delete \
    && rm -f /usr/lib/x86_64-linux-gnu/libcudnn_static_v*.a \ 
    && ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/local/cuda/lib64/stubs/libcuda.so.1 \
    && echo "/usr/local/cuda/lib64/stubs" > /etc/ld.so.conf.d/z-cuda-stubs.conf \
    && ldconfig

# Install miniconda
RUN MINICONDA="Miniconda3-latest-Linux-x86_64.sh" \
    && wget --quiet https://repo.continuum.io/miniconda/$MINICONDA \
    && bash $MINICONDA -b -p /miniconda \
    && rm -f $MINICONDA

ENV PATH /miniconda/bin:$PATH

COPY environment.yml .
RUN conda env create -f environment.yml
