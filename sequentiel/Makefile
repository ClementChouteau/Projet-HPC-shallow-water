SRC:=export.cpp forward.cpp init.cpp memory.cpp parse_args.cpp shalw.cpp
OBJ:=$(addprefix obj/, $(SRC:.cpp=.o))
SRC:=$(addprefix src/, $(SRC))
INCLUDE:=include
TARGET=shalw
ARGS=
EXPORT_PATH=/tmp
SIZE_TEST:=512
T_TEST:=40
TEST_NAME:=shalw_$(SIZE_TEST)x$(SIZE_TEST)_T$(T_TEST)

CXX=mpicxx
CXXFLAGS=-Wall -Wextra -Werror -O3 -I $(INCLUDE)
LDFLAGS=-lboost_program_options -lm -O3

all: rep bin/$(TARGET)

rep:
	@if [ ! -d "bin" ]; then mkdir bin; fi
	@if [ ! -d "obj" ]; then mkdir obj; fi

bin/$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) $(ARGS)

exec: rep bin/$(TARGET)
	bin/$(TARGET) $(ARGS)

obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

test1: rep bin/$(TARGET)
	bin/$(TARGET) --export --export-path $(EXPORT_PATH)

test2: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 8192 -y 8192 -t 20

testSave: all
	bin/$(TARGET) -x $(SIZE_TEST) -y $(SIZE_TEST) -t $(T_TEST) --export --export-path ./bin/
	../save_last.py ./bin/$(TEST_NAME).sav

clean:
	rm -f obj/*.o

mrproper: clean
	rm -f bin/*

.PHONY: rep exec test1 test2 testSave clean mrproper
