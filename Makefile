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
	rm -f obj/*.o

mrproper : clean
	rm -f bin/shalw_seq
	rm -f bin/shalw_bandes
	rm -f bin/shalw_bandes_nblock
