######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)


CFLAGS= -Wall -std=c++0x -O3 -DVERSION=\"$(GIT_VERSION)\"



all: pi

pi:
	$(CC) $(CFLAGS) -I fastahack fastahack/*o src/pi.cpp -o PI