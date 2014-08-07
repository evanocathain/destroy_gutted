CC = gcc
CFLAGS = -c -O3 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64
INCL = -I./ #-I/home/pulsar/pgplot/
LIBS = -lm
#LIBS = -L/usr/X11/lib -lX11 -L/home/pulsar/pgplot/ -lg2c -L/usr/lib/ -lpng -L/usr/lib64/ -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/ -lm

all: destroy strip clean

destroy: destroy.o get_args.o read_data.o strings_equal.o sp_search.o
	$(CC) destroy.o get_args.o read_data.o strings_equal.o sp_search.o -o destroy $(LIBS)

destroy.o: destroy.c get_args.c read_data.c strings_equal.c
	$(CC) $(CFLAGS) $(INCL) destroy.c

get_args.o: get_args.c
	$(CC) $(CFLAGS) $(INCL) get_args.c

read_data.o: read_data.c
	$(CC) $(CFLAGS) $(INCL) read_data.c

strings_equal.o: strings_equal.c
	$(CC) $(CFLAGS) $(INCL) strings_equal.c

sp_search.o: sp_search.c
	$(CC) $(CFLAGS) $(INCL) sp_search.c

strip:
	strip destroy

clean:
	rm -f *~ *.o
