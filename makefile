#
#   arraygen makefile
#

CC=gcc -lm -g -Wall -O3


all: 
	$(CC) arraygen.c -o arraygen

clean: 
	rm arraygen 
