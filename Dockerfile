FROM ubuntu:18.04
COPY ./BOLT /workspace/BOLT
WORKDIR /workspace/

RUN apt-get update
RUN apt-get install -y g++
RUN apt-get install -y libblas-dev liblapack-dev

RUN g++ --std=c++0x BOLT/src/* -llapack -o BOLT/run
