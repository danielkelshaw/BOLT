FROM ubuntu:18.04
COPY ./BOLT /workspace/BOLT
WORKDIR /workspace/

RUN apt-get update
RUN apt-get install -y g++
RUN apt-get install -y libblas-dev liblapack-dev
RUN apt-get install -y make

RUN make --directory /workspace/BOLT/
