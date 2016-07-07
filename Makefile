CC = g++-6
CFLAGS = -g -Wall -O3 -std=c++14 -fPIC -Ic:/gsl -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgslcblas 
LDFLAGS= -Lc:/gsl
LIBS= -lgsl
SRCS = plot.cpp gaussfit.cpp kde.cpp kgpmain.cpp
OBJS = $(SRCS:.cpp=.o)

all: kgp

$(OBJS): %.o : %.h

kgpmain.o: gaussfit.h plot.h kde.h kgpmain.cpp

.c.o:
	$(CC) -c $< $(CFLAGS)

kgp: $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(CFLAGS)

lib: $(OBJS)
	$(CC) -shared -o kgp.so $^  $(LDFLAGS) $(LIBS) $(CFLAGS)

clean:
	rm -f $(OBJS) kgp kgp.so