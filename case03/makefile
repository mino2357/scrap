CXX       = g++
CXXFLAGS  = -std=c++17 -O2 -Wall -Wextra
LDFLAGS   =
TARGET    = sim
SOURCES   = main.cpp
HEADERS   = porous_reactor.hpp
PLOT_FILE = plot_live.gp

.PHONY: all clean live run plot

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@

# スナップショット生成（最後に conc.dat / temp.dat が残る）
run: $(TARGET)
	./$(TARGET)

# ライブ可視化（ファイル監視）
live: $(TARGET) plot_live.gp
	@trap 'kill -TERM $$(cat sim.pid) 2>/dev/null || true; rm -f sim.pid' INT TERM EXIT; \
	./$(TARGET) & echo $$! > sim.pid; \
	gnuplot -persist plot_live.gp

clean:
	rm -f $(TARGET) conc.dat temp.dat conc.tmp temp.tmp sim.pid
