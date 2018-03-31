SRC:=export.cpp forward.cpp init.cpp memory.cpp parse_args.cpp shalw.cpp
OBJ:=$(addprefix obj/, $(SRC:.cpp=.o))
SRC:=$(addprefix src/, $(SRC))
INCLUDE:=include
TARGET=shalw
ARGS=
EXPORT_PATH=/tmp

CXX=mpicxx
CXXFLAGS=-Iinc -Wall -O3 -I $(INCLUDE)
LDFLAGS=-lboost_program_options -lm -O3

all: rep bin/$(TARGET)

rep:
	@if [ ! -d "bin" ]; then mkdir bin; fi
	@if [ ! -d "obj" ]; then mkdir obj; fi

bin/$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) $(ARGS)

exec: all
	bin/$(TARGET) $(ARGS)

obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

test1: all
	bin/$(TARGET) --export --export-path $(EXPORT_PATH)

testSave: all
	bin/$(TARGET) -x 512 -y 512 -t 40 --export --export-path ./bin/

clean:
	rm -f obj/*.o

mrproper: clean
	rm -f bin/*