CC	= g++

CFLAGS	= -O3 -o

SEQPROG	= rainfall_seq

PARAPROG = rainfall_pt

PROG = $(SEQPROG) $(PARAPROG)

LIB = -lpthread

all: $(PROG)

$(SEQPROG): rainfall_seq.cpp
	$(CC) $(CFLAGS) $(SEQPROG) rainfall_seq.cpp
$(PARAPROG):rainfall_pt.cpp
	$(CC) $(CFLAGS) $(PARAPROG) rainfall_pt.cpp $(LIB)
clean:
	rm -f *.o $(PROG)
