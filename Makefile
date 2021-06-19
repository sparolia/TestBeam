
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I../main/include/
#-I/home/pitzl/GBL/V01-17-00/cpp/include

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pitzl/GBL/V01-17-00/cpp/lib/

scope53m: scope53m.cc
	g++ $(CXXFLAGS) scope53m.cc -o scope53m \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: scope53m'

edg53: edg53.cc
	g++ $(CXXFLAGS) edg53.cc -o edg53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: edg53'

scope53: scope53.cc
	g++ $(CXXFLAGS) scope53.cc -o scope53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: scope53'

tele: tele.cc
	g++ tele.cc $(CXXFLAGS) -fopenmp -o tele \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: tele'

ed53: ed53.cc
	g++ $(CXXFLAGS) ed53.cc -o ed53 \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: ed53'

ed3d: ed3d.cc
	g++ $(CXXFLAGS) ed3d.cc -o ed3d \
	$(ROOTLIBS) -L../lib -lEUDAQ
	@echo 'done: ed3d'
