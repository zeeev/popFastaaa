######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)


CFLAGS= -Wall -std=c++0x -O3 -DVERSION=\"$(GIT_VERSION)\"


all: fastaH pi snp indel

fastaH:
	cd fastahack/ && make

pi:
	$(CC) $(CFLAGS) -I fastahack fastahack/*o src/pi.cpp -o PI
snp:
	$(CC) $(CFLAGS) -I fastahack fastahack/*o src/snp.cpp -o SNP
indel:
	$(CC) $(CFLAGS) -I fastahack fastahack/*o src/indel.cpp -o INDEL