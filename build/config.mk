# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))

PREFIX=/home/aschneider/programs/GOLEMSPACE/local

# General Settings
VERSION=1.0.0
CXXFLAGS= -std=c++11 -O2 -g -fPIC

# Libraries
CXXFLAGS+=-I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/include -I/usr/include
LDFLAGS+=-fPIC -L/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib -lboost_system  -lboost_iostreams -L/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib -lhdf5_hl -lhdf5 -L/usr/lib -lrt -lz -ldl -lm -Wl,-rpath -Wl,/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib


