# Container image that runs your code
FROM jcsda/docker-gnu-openmpi-dev:latest AS base
RUN apt-get update \
    && apt-get -y install libhdf5-serial-dev \
    && apt-get -y install libnetcdf-dev libnetcdff-dev libnetcdf-c++4-dev \
    && apt-get -y install lz4

# Something like the below pattern would be quicker than using github-actions
#caching, however I would likely need to slim down the jcsda image as it is
#rather fat. I need a different version of eckit that contains LZ4, but I don't
#think I need any of the python stuff for example.

#FROM base AS eckitimage
#RUN cd /root \
#    && git clone https://github.com/ecmwf/eckit.git \
#    && mkdir eckit_build \
#    && cd eckit_build \
#    && /usr/local/bin/ecbuild --prefix=/usr/local -- ../eckit \
#    && make \
#    && make install \
#    && cd ../ \
#    && rm -rf eckit_build

