
export CXX="g++"
export CXXFLAGS="-std=c++11 -O0 -g -I../ -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/include -I/usr/include"
export LDFLAGS="-L../build -lPhysTools -L/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib -lboost_system  -lboost_iostreams -L/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib -lhdf5_hl -lhdf5 -L/usr/lib -lrt -lz -ldl -lm -Wl,-rpath -Wl,/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib"

export LD_LIBRARY_PATH="../lib:/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_6_x86_64/lib"
