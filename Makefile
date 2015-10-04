JC = javac
JFLAGS = -Xlint:all

SRC = $(wildcard *.java)

EXE = ducks
CLASSES = $(wildcard *.java)
FIFO = player2server

all: $(EXE)

$(EXE): $(SRC)
	$(JC) $(JFLAGS) $^ 

game:
	rm -f $(FIFO) && mkfifo $(FIFO)
	java Main server < $(FIFO) | java Main verbose > $(FIFO)

clean:
	rm -f $(EXE)
	rm -f $(CLASSES)
	rm -f $(FIFO)
