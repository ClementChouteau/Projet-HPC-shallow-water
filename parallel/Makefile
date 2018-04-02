SRC:=export.cpp forward.cpp init.cpp memory.cpp parse_args.cpp shalw.cpp
OBJ:=$(addprefix obj/, $(SRC:.cpp=.o))
SRC:=$(addprefix src/, $(SRC))
INCLUDE:=include
TARGET:=shalw
ARGS:=
NODES:=2
HOSTFILE:=../hostfile

CXX=mpicxx
CXXFLAGS=-Iinc -Wall -O3 -I $(INCLUDE) -g
LDFLAGS=-lboost_program_options -lm
EX=mpirun
EXFLAGS=-n $(NODES) -hostfile $(HOSTFILE) --map-by node

all: rep bin/$(TARGET)

rep:
	@if [ ! -d "bin" ]; then mkdir bin; fi
	@if [ ! -d "obj" ]; then mkdir obj; fi

bin/$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

exec: all
	$(EX) $(EXFLAGS) bin/$(TARGET) $(ARGS)

obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

test2: all
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 8192 -y 8192 -t 20

test3: all
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 512 -y 512 -t 40 --export --export-path ~/

testSave: all
	$(EX) -n 4 -hostfile $(HOSTFILE) --map-by node bin/$(TARGET) -x 512 -y 512 -t 40 --block --export --export-path ./bin/
	../save_last.py ./bin/shalw_512x512_T40.sav

io: all
	$(EX) -n 4 -hostfile $(HOSTFILE) --map-by node bin/$(TARGET) -x 512 -y 512 -t 40 --block --export --export-path ./bin/
	../save_last.py bin/shalw_512x512_T40.sav
	gimp ./bin/shalw_512x512_T40.png

clean:
	rm -f obj/*.o

mrproper: clean
	rm -f bin/*