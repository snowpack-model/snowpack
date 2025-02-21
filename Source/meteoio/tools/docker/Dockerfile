FROM alpine AS builder
RUN apk update && apk add  --no-cache git g++ make cmake netcdf-dev proj-dev curl
RUN mkdir -p /root/usr && \
    mkdir -p /root/src && \
    cd /root/src && \
    git clone https://gitlabext.wsl.ch/snow-models/meteoio.git && \
    cd meteoio && \
    cmake -D CMAKE_INSTALL_PREFIX:PATH=/root/usr -D DEST:STRING=optimized -D PLUGIN_NETCDFIO:BOOL=ON -D PLUGIN_DBO:BOOL=ON -D VERSION_FROM_GIT:BOOL=ON -D PROJ:BOOL=ON && \
    make -j 2 && \
    make install

FROM alpine AS meteoio
ENV PATH=/root/usr/bin:$PATH
ENV METEOIO_VERSION=2.10.0
ENV RELEASE_DATE=2021-10-05
LABEL slf.ch.version="${METEOIO_VERSION}" com.example.release-date="${RELEASE_DATE}" vendor="slf.ch" licence="LGPL-v3"
RUN apk update && apk add  --no-cache netcdf proj curl
COPY --from=builder /root/usr /root/usr
CMD ["bin/sh"]
