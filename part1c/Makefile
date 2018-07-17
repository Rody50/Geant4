# project name (generates executable with this name)
TARGET = Diagonalization

CC=g++
# compiling flags
CFLAGS=-c -std=c++11 -O3

# linking flags
LDFLAGS=-larmadillo

# libraries
LIBS=-lm -lgsl -lgslcblas -larmadillo -llapack

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