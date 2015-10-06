#!/bin/sh

mpicxx -std=c++11 -o mpi-wrapper mpi-wrapper.cpp
mv mpi-wrapper ~/bin/
