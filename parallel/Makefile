SRC:=export.cpp forward.cpp init.cpp memory.cpp parse_args.cpp shalw.cpp
OBJ:=$(addprefix obj/, $(SRC:.cpp=.o))
SRC:=$(addprefix src/, $(SRC))
INCLUDE:=include
TARGET:=shalw
ARGS:=
MODE:=--async --block
NODES:=4
HOSTFILE:=../hostfile
SIZE_TEST:=512
T_TEST:=40
EXPORT_STEP:=1
TEST_NAME:=shalw_$(SIZE_TEST)x$(SIZE_TEST)_T$(T_TEST)

CXX=mpicxx
CXXFLAGS=-Wall -Werror -Wextra -O3 -I $(INCLUDE) -fopenmp
LDFLAGS=-lboost_program_options -lm -fopenmp
EX=mpirun
EXFLAGS=-n $(NODES) -hostfile $(HOSTFILE) --map-by node

all: rep bin/$(TARGET)

rep:
	@if [ ! -d "bin" ]; then mkdir bin; fi
	@if [ ! -d "obj" ]; then mkdir obj; fi

bin/$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

exec: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) $(ARGS) $(MODE)

obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

test2: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 8192 -y 8192 -t 20 $(MODE)

test3: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 512 -y 512 -t 40 --export --export-path ~/ --export-step $(EXPORT_STEP)

big_test: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 32768 -y 32768 -t 40 $(MODE)

small_test: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 2048 -y 2048 -t 20 $(MODE)

testSave: rep bin/$(TARGET)
	-rm ./bin/$(TEST_NAME).sav	
	$(EX) -n 4 -hostfile $(HOSTFILE) --map-by node bin/$(TARGET) -x $(SIZE_TEST) -y $(SIZE_TEST) -t $(T_TEST) $(MODE) --export --export-path ./bin/ --export-step $(EXPORT_STEP)
	../save_last.py ./bin/$(TEST_NAME).sav

io: rep bin/$(TARGET)
	-rm ./bin/$(TEST_NAME).png
	-rm ./bin/$(TEST_NAME).sav
	$(EX) -n 4 -hostfile $(HOSTFILE) --map-by node bin/$(TARGET) -x $(SIZE_TEST) -y $(SIZE_TEST) -t $(T_TEST) $(MODE) --export --export-path ./bin/ --export-step $(EXPORT_STEP)
	../save_last.py bin/$(TEST_NAME).sav
	gimp ./bin/$(TEST_NAME).png

test_hybride: rep bin/$(TARGET)
	$(EX) $(EXFLAGS) bin/$(TARGET) -x 2048 -y 2048 -t 20 --async --block --hybride

clean:
	rm -f obj/*.o

mrproper: clean
	rm -f bin/*

.PHONY: rep exec test2 test3 testSave io clean mrproper
