CC = g++
CXXFLAGS = -Wall -std=c++0x 

EXE = ducks
OBJFILES = main.o hmm.o Player.o GameServer.o Client.o

all: $(EXE)

$(EXE): $(OBJFILES)
	mkfifo player2server
	$(LINK.cpp) $(LOADLIBES) $(LDLIBS) $^ -o $@

gameself:
	mkfifo player2server
	./$(EXE) server < player2server | ./$(EXE) verbose > player2server

clean:
	rm -f player2server && mkfifo player2server
	rm -f $(EXE)
	rm -f $(OBJFILES)
	rm -f pipe
