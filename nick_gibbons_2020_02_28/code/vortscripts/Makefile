# for Mac OS
# make all "CFLAGS=-Wall -DMACOSX"

CC= gcc

OPT = -O3

ifeq ($(MAKECMDGOALS), debug)
    OPT = -O0 -g
endif

FLAGS = -Wall -fPIC -I.
CFLAGS = ${FLAGS} ${OPT}

LIBS= -lm -llapack

# These seem to hapen automatically somehow. 
#nmsimplex.o: nmsimplex.c nmsimplex.h
#	${CC} ${CFLAGS} ${OPT} nmsimplex.c ${LIBS}
#
#vortid.o: vortid.c
#	${CC} ${CFLAGS} ${OPT} vortid.c ${LIBS}

debug: vortid.o nmsimplex.o
	${CC} ${CFLAGS} vortid.o nmsimplex.o ${LIBS} -o test 

shared: vortid.o nmsimplex.o
	${CC} ${CFLAGS} -shared  vortid.o nmsimplex.o ${LIBS} -o libvortid.so


clean:
	rm -f ./*.o 
