CC=gcc
FLAG= -lgmp -Wall -W -g
OBJ= main.o algo.o
JOB=./include/*.h
SRC=./src/

all:fact clean doc

fact: ${OBJ}
	${CC} -o ./fact ${OBJ} ${FLAG}

algo.o: ${SRC}algo.c ${JOB}
	${CC} -c ${SRC}algo.c ${FLAG}

main.o: ${SRC}main.c ${JOB}
	${CC} -c ${SRC}main.c ${JOB} ${FLAG}

clean:
	rm -rf *.o
	rm -rf ./include/*.h.gch
	rm -rf doc
doc:
	doxygen
