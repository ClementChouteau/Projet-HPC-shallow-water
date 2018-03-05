CXX=mpicxx -Wall -O3

CPPFLAGS=-Iinc

LDFLAGS=-lboost_program_options -lm

FILES=export forward init memory parse_args shalw
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/shalw_seq

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LDFLAGS)

obj/%.o : src_seq/%.cpp
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

obj/%.o : src_seq/%.c
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

clean :
	rm -f bin/* obj/*
