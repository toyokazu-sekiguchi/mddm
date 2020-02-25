CC = gcc
AR = ar rvs
#CCFLAG = -O2 -DCAMB -DHYRECPATH=\"$(PWD)/\" 
CCFLAG = -Wall -fPIC -DHYRECPATH=\"$(PWD)/\" 
LDFLAG = -O2

INC_PY = -I/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m
LIB_PY = -L/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib -lpython3.7m
#-I/anaconda3/include/python3.6m
#-L/anaconda3/lib -lpython3.6m
#-Wunused-command-line-argument
#-I/anaconda3/include -arch x86_64

%.o: %.c 
	$(CC) $(CCFLAG) -c $*.c -o $*.o

HYREC_SRC = hyrectools.o helium.o hydrogen.o history.o 
HYREC_EXE = hyrec.o 

HYREC_SRC2 = hyrectools.c helium.c hydrogen.c history.c hyrec.c

default: libhyrec.a

clean: 
	rm *.o

hyrec: $(HYREC_SRC2) 
	$(CC) $(LDFLAG) -o hyrec  $(HYREC_SRC2) -lm

pyrec: $(HYREC_SRC) 
	$(CC) $(CCFLAG) -c pyrec.c -o pyrec.o
	swig -python pyrec.i
	$(CC) $(CCFLAG) $(INC_PY) -c pyrec_wrap.c -o pyrec_wrap.o
	$(CC) $(LIB_PY) -shared $(HYREC_SRC) pyrec.o pyrec_wrap.o -o _pyrec.so

libhyrec.a: $(HYREC_SRC)
	$(AR) $@ $(HYREC_SRC)
