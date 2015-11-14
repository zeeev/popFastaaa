######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
CFLAGS= -Wall -std=c++0x -O3

all: pi

pi:
	$(CC) $(CFLAGS) -I fastahack fastahack/*o src/pi.cpp -o PI