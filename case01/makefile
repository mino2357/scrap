CXX = /usr/bin/g++
TARGET = main
CXXFLAGS = -Wall -Wextra -O2 -march=native -std=c++2b 
LDLFLAGS = -lstdc++
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS)

run: all
	./$(TARGET)

clean:
	rm $(TARGET)
