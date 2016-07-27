# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

all: HoughFilter

HoughFilter: HoughFilter.o Settings.o Stub.o Cell.o Hough.o 
	$(CC) HoughFilter.o Settings.o Stub.o Hough.o Cell.o $(ROOTFLAGS) $(ROOTLIBS) -L/usr/local/lib -lboost_filesystem -o HoughFilter

HoughFilter.o: src/HoughFilter.cc
	$(CC) $(CFLAGS) $(ROOTFLAGS) src/HoughFilter.cc
Settings.o: src/Settings.cc 
	$(CC) $(CFLAGS) $(ROOTFLAGS) src/Settings.cc
Stub.o: src/Stub.cc
	$(CC) $(CFLAGS) $(ROOTFLAGS) src/Stub.cc
Hough.o: src/Hough.cc 
	$(CC) $(CFLAGS) $(ROOTFLAGS) src/Hough.cc
Cell.o: src/Cell.cc 
	$(CC) $(CFLAGS) $(ROOTFLAGS) src/Cell.cc

clean:
	rm *o HoughFilter
