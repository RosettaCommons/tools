#
# Makefile for C version of MaxSprout
#
# Chris Dodge - EBI - Oct 96
# Ville Silventoinen - EBI - May 2002 (adapted to Linux)
#
CC = gcc
#CC = pgcc
CFLAGS = -g -O
#CFLAGS = -O3 -pipe -mcpu=pentiumpro -march=pentiumpro -fomit-frame-pointer
FFLAGS = -C -O

all: readbrk buildbackbone torso2

buildbackbone: buildbackbone.o u3b.o
	${CC} $(CFLAGS) -o buildbackbone buildbackbone.o u3b.o -lg2c -lm -lp2c

readbrk: readbrk.o u3b.o
	${CC} $(CFLAGS) -o readbrk readbrk.o u3b.o -lg2c -lm -lp2c

readbrk.o: readbrk.c
	${CC} $(CFLAGS) -c readbrk.c

buildbackbone.o: buildbackbone.c
	${CC} $(CFLAGS) -c buildbackbone.c

u3b.o:	u3b.f
	f77 $(FFLAGS) -c u3b.f

torso2: ran.o torso.o torsotool.o u3b.o
	f77 $(FFLAGS) -o torso2 ran.o torso.o torsotool.o u3b.o

clean:
	rm -f *.o buildbackbone readbrk torso2
