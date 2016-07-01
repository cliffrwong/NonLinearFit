CC = g++-5
CFLAGS = -g -Wall -O3 -fPIC -Ic:/gsl -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgslcblas -lsndfile
LDFLAGS= -Lc:/gsl
LIBS= -lgsl
SRCS = wavproc.c gaussfit.c
OBJS = $(SRCS:.cpp=.o)

all: wavproc

$(OBJS): %.o : %.h

wavproc.o: gaussfit.hpp

.c.o:
	$(CC) -c $< $(CFLAGS)

wavproc: $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(CFLAGS)

lib: $(OBJS) wavproc.o
	$(CC) -shared -o wavproc.so $^  $(LDFLAGS) $(LIBS) $(CFLAGS)

clean:
	rm -f $(OBJS) wavproc wavproc.so