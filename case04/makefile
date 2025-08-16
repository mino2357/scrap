CXX=g++
CXXFLAGS=-std=c++17 -O2 -Wall -Wextra
TARGET=chem
SOURCES=main.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
