# Container image that runs your code
FROM jcsda/docker-gnu-openmpi-dev:revert AS base

RUN apt-get update \
    && apt-get -y install libhdf5-serial-dev \
    && apt-get -y install libnetcdf-dev libnetcdff-dev libnetcdf-c++4-dev \
    && cd /root \
    && git clone https://github.com/lz4/lz4.git \
    && cd lz4 \
    && git checkout tags/v1.9.3 \
    && make \
    && make --eval prefix=/opt install \
    && export LZ4_ROOT=/opt

# latest JCSDA images use spack to install the environment
#source /etc/profile.d/z10_spack_environment.sh

# Something like the below pattern would be quicker than using github-actions
#caching, however I would likely need to slim down the jcsda image as it is
#rather fat. I need a different version of eckit that contains LZ4, but I don't
#think I need any of the python stuff for example.

# NOTE: to use this I will need to remove eckit/fckit from CMakeLists and ensure
#       library is available at ecbuild time 
#FROM base AS eckitimage
#RUN cd /root \
#  && git clone https://github.com/ecmwf/eckit.git \
#  && cd eckit && git checkout tags/1.19.0 && cd ../ \
#  && mkdir eckit_build \
#  && cd eckit_build \
#  && /usr/local/bin/ecbuild --prefix=/usr/local -- ../eckit \
#  && make \
#  && make install \
#  && cd ../ \
#  && rm -rf eckit_build

#FROM base AS fckitimage
#RUN cd /root \
#  && git clone https://github.com/ecmwf/fckit.git \
#  && cd fckit && git checkout tags/0.9.5 && cd ../ \
#  && mkdir fckit_build \
#  && cd fckit_build \
#  && /usr/local/bin/ecbuild --prefix=/usr/local -- ../fckit \
#  && make \
#  && make install \
#  && cd ../ \
#  && rm -rf fckit_build
