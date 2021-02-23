CC=mpicc -std=gnu99 -Wall -Wno-unused-variable -O3 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

AR=ar

INCLUDE=../include/*.h
INCLDDIR= -I../include/ -I./ -I/mnt/c/Users/夏/Desktop/lime-1.3.2/include/ -I/mnt/c/Users/夏/Desktop/DDalphaAMG/include/

TARGET=./libqcd.a

all: $(TARGET)

SOURCES=\
	qcd_blas\
	qcd_communication\
	qcd_gamma\
	qcd_init\
	qcd_io\
	qcd_observables\
	qcd_smearing\
	qcd_smearing_special\
	qcd_blas_special\
        qcd_func


$(TARGET):		${addsuffix .o, $(SOURCES)} $(INCLUDE)
			$(AR) -r $@ ${addsuffix .o, $(SOURCES)}
			ranlib $@

%.o:			%.c $(INCLUDE)
			$(CC) -c $< $(INCLDDIR)
clean:		
		rm -f *.o libqcd.a
