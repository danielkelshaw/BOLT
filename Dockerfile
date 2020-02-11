FROM ubuntu:18.04
COPY ./BOLT /workspace/BOLT
WORKDIR /workspace/

RUN mkdir /workspace/BOLT/obj

RUN apt-get update
RUN apt-get install -y sudo

RUN sudo apt-get install -y g++
RUN sudo apt-get install -y make

RUN sudo apt-get install -y libblas-dev 
RUN sudo apt-get install -y liblapack-dev
RUN sudo apt-get install -y libboost-all-dev

RUN make --directory BOLT/
