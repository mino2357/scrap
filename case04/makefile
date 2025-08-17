.RECIPEPREFIX := >
CXX=g++
CXXFLAGS=-std=c++17 -O2 -Wall -Wextra
TARGET=chem
SOURCES=main.cpp

TEST_TARGET=test_harmonic
TEST_SOURCES=test_harmonic.cpp

.PHONY: all clean run test

all: $(TARGET)

$(TARGET): $(SOURCES)
>$(CXX) $(CXXFLAGS) $(SOURCES) -o $@

$(TEST_TARGET): $(TEST_SOURCES)
>$(CXX) $(CXXFLAGS) $(TEST_SOURCES) -o $@

run: $(TARGET)
	python3 plot_thermo.py
	./$(TARGET) rk78
	gnuplot plot.gp

clean:
>rm -f $(TARGET) $(TEST_TARGET)
