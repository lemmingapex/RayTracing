CFLAGS = -O3 -funroll-loops -g -Wall -ansi
#CFLAGS = -O1 -g -Wall -ansi
#CFLAGS = -g -Wall -ansi
CC = g++
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)
PROG = RayTrace

all: $(SRCS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@ 

depend:
	makedepend -Y $(SRCS)

clean:
	rm $(OBJS) $(PROG)
