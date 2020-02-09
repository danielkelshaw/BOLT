# BOLT makefile

CC=g++
CFLAGS=--std=c++0x -Wall

EXE=BOLT

DIR=.
SOURCEDIR=$(DIR)/src
HEADERDIR=$(DIR)/inc
OBJECTDIR=$(DIR)/obj

SRCS:=$(wildcard $(SOURCEDIR)/*.cpp)
OBJS:=$(addprefix $(OBJECTDIR)/, $(notdir $(SRCS:.cpp=.o)))

LIB=-llapack -lboost_system -lboost_filesystem

# Build BOLT
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

# Build object files
$(OBJS): $(OBJECTDIR)/%.o : $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf $(EXE) $(OBJECTDIR) Results && mkdir $(OBJECTDIR) && mkdir Results && mkdir Results/VTK

.PHONY: depend
depend:
	makedepend -Y -- $(CFLAGS) -- $(SRCS)

