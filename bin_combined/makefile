CC=mpicc

CFLAGS=-I../include -std=gnu99 -O3 -I/mnt/c/Users/夏/Desktop/DDalphaAMG/include/
LDFLAGS=-L../lib -L/mnt/c/Users/夏/Desktop/lime-1.3.2/lib  -L/mnt/c/Users/夏/Desktop/DDalphaAMG/lib/  #-L/mnt/c/Users/夏/Desktop/package/DDalphaAMG-master/build/  # -L/PATH_TO_LIME/
LIBS=-lqcd -llime -lm -lDDalphaAMG -fopenmp

.PHONY: clean\
	cleanall\
	lib

#TARGETS=show_conf

TARGETS=final\
       # twopdelta  \
	threeponeplusgamma\
	sourceseqoneplusgamma\
	sourceseqdnoneplusgamma\
	sourceseqoneforvector\
	sourceseqdnoneforvector\
	threeponeforvector\
        sourceseqdnoneplusgammafixcur\
        sourcemom\
        quaipdf	\
	sourcemomshape

all: lib ${addsuffix .exe, $(TARGETS)}

lib: ../lib/libqcd.a
../lib/libqcd.a:
	cd ../lib &&\
	make

%.exe: %.o ../lib/libqcd.a
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

%.o: %.c 
	$(CC) $(CFLAGS) -c $<

clean:		
	rm -vf ${addsuffix .o, $(TARGETS)}

cleanall: clean
	rm -vf ${addsuffix .exe, $(TARGETS)} &&\
	cd ../lib/ &&\
	make clean
