JC = javac
JFLAGS = -Xlint:all

SRC = $(wildcard *.java)

EXE = ducks
CLASSES = $(wildcard *.class)
FIFO = player2server

INPUTFILE = SouthEmissions

exe: $(EXE)

$(EXE): $(SRC)
	$(JC) $(JFLAGS) $^ 

game: fifo
	java Main server < $(FIFO) | java Main verbose > $(FIFO)

north: fifo
	java Main load NorthEmissions.in server < $(FIFO) | java Main load NorthEmissionsOpponent.in verbose > $(FIFO)

south: fifo
	java Main load SouthEmissions.in server < $(FIFO) | java Main load SouthEmissionsOpponent.in verbose > $(FIFO)

west: fifo
	java Main load WestEmissions.in server < $(FIFO) | java Main load WestEmissionsOpponent.in verbose > $(FIFO)

east: fifo
	java Main load EastEmissions.in server < $(FIFO) | java Main load EastEmissionsOpponent.in verbose > $(FIFO)

paradise: fifo
	java Main load ParadiseEmissions.in server < $(FIFO) | java Main load ParadiseEmissionsOpponent.in verbose > $(FIFO)

all: north south west east paradise

fifo:
	rm -f $(FIFO) && mkfifo $(FIFO)

clean:
	rm -f $(EXE)
	rm -f $(CLASSES)
	rm -f $(FIFO)
