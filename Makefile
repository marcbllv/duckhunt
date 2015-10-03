CC = g++
CXXFLAGS = -Wall -std=c++0x 

EXE = ducks
OBJFILES = main.o hmm.o Player.o GameServer.o Client.o
FIFO = player2server

all: $(EXE)

$(EXE): $(OBJFILES)
	$(LINK.cpp) $(LOADLIBES) $(LDLIBS) $^ -o $@

game:
	rm -f $(FIFO) && mkfifo $(FIFO)
	./$(EXE) server < $(FIFO) | ./$(EXE) verbose > $(FIFO)

clean:
	rm -f $(EXE)
	rm -f $(OBJFILES)
	rm -f $(FIFO)
