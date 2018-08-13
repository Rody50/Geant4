# project name (generates executable with this name)
TARGET = talent

CC=g++
# compiling flags
CFLAGS=-c -std=c++11 -O3 -g -ggdb

# linking flags
LDFLAGS=-larmadillo

# libraries
LIBS=-lm -larmadillo -llapack -lblas

# compilation options for ROOT
#RO=$(shell root-config --cflags) $(shell root-config --libs)
#ROOTINC=$(shell root-config --incdir)
#LIBS+=$(RO)

# specific directories for each extension
SRCDIR=./src
OBJDIR=./src/obj
BINDIR=./bin

SOURCES=$(wildcard $(SRCDIR)/*.cpp) 
HEADERS=$(wildcard $(SRCDIR)/*.h)
OBJECTS=$(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(TARGET)
