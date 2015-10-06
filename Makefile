JC = javac
JFLAGS = -Xlint:all

SRC = $(wildcard *.java)

EXE = ducks
CLASSES = $(wildcard *.class)
FIFO = player2server

INPUTFILE = SouthEmissions

all: $(EXE)

$(EXE): $(SRC)
	$(JC) $(JFLAGS) $^ 

game:
	rm -f $(FIFO) && mkfifo $(FIFO)
	java Main load $(INPUTFILE).in server < $(FIFO) | java Main load $(INPUTFILE)Opponent.in verbose > $(FIFO)

clean:
	rm -f $(EXE)
	rm -f $(CLASSES)
	rm -f $(FIFO)
