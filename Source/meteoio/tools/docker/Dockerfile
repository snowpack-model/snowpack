FROM alpine AS builder
RUN apk add  --no-cache git g++ make cmake netcdf-dev proj-dev
RUN mkdir -p /root/usr && \
    mkdir -p /root/src && \
    cd /root/src && \
    git clone https://gitlabext.wsl.ch/snow-models/meteoio.git && \
    cd meteoio && \
    cmake -D CMAKE_INSTALL_PREFIX:PATH=/root/usr -D DEST:STRING=optimized -D PLUGIN_NETCDFIO:BOOL=ON -D VERSION_FROM_GIT:BOOL=ON -D PROJ:BOOL=ON && \
    make && \
    make install

FROM alpine
LABEL slf.ch.version="2.9.0" com.example.release-date="2021-10-05" vendor="slf.ch" licence="LGPL-v3"
RUN apk add  --no-cache netcdf proj
RUN mkdir -p /root/usr && \
    cd /root/usr
COPY --from=builder /root/usr /root/usr
ENV PATH=/root/usr/bin:$PATH
