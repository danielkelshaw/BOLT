FROM ubuntu:18.04

RUN mkdir -p /workspace/

COPY ./BOLT /workspace/BOLT
RUN mkdir -p /workspace/BOLT/obj

WORKDIR /workspace/

RUN apt-get update && apt-get install -y \
	g++ \
	make \
	libblas-dev \
	liblapack-dev \
	libboost-all-dev \

RUN make --directory BOLT/

