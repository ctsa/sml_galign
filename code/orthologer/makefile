
INC_DIR   = 
LIB_DIR   = 

OBJS	= orthologer.o 

# General flags
CXX      =  g++
MFLAGS   = -O3
DBFLAGS  = -ggdb -D_GLIBCXX_DEBUG
WFLAGS   = -Wall
CXXFLAGS = $(MFLAGS) $(WFLAGS)
LIBS	 =
PROG     = orthologer

all:	$(OBJS)
	$(CXX) $(CXXFLAGS) -static $(OBJS) $(LIBS) -o $(PROG) 
	-rm -f core *.o

clean:	
	-rm -f core *.o

