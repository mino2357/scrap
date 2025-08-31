.RECIPEPREFIX := >
CXX=g++
CXXFLAGS=-std=c++17 -Ofast -march=native -funroll-loops -fomit-frame-pointer -Wall -Wextra
LDLIBS=-lquadmath
TARGET=three_body
SOURCES=main.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
>$(CXX) $(CXXFLAGS) $(SOURCES) $(LDLIBS) -o $@

run: $(TARGET)
>./$(TARGET)

clean:
>rm -f $(TARGET) case05.dat
