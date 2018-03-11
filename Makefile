CXX=mpicxx

CPPFLAGS=-Iinc -Wall -O3

LDFLAGS=-lboost_program_options -lm

FILES=export forward init memory parse_args shalw
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/shalw_bandes_nblock

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LDFLAGS)

obj/%.o : src_bandes_nblock/%.cpp
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

obj/%.o : src_bandes_nblock/%.c
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

clean :
	rm -f obj/*.o

test : all
	./check_last.sh


mrproper : clean
	rm -f bin/shalw_seq
	rm -f bin/shalw_bandes
	rm -f bin/shalw_bandes_nblock
