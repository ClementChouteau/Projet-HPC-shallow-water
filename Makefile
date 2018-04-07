VERSIONS:=sequentiel parallel
NODES:=16
ROOM:=401
SLOTS:=4
TEST:=test2
TIME:=time -f "%e %C" -ao times.csv

all: sequentiel parallel

parallel: hostfile
	$(MAKE) -C parallel

sequentiel: hostfile
	$(MAKE) -C sequentiel

times: hostfile parallel
	@echo -n "" > times.csv
	@$(TIME) $(MAKE) NODES=$(NODES) -C parallel $(TEST)
	@$(TIME) $(MAKE) NODES=$(NODES) MODE=--async -C parallel $(TEST)
	@$(TIME) $(MAKE) NODES=$(NODES) MODE=--block -C parallel $(TEST)
	@$(TIME) $(MAKE) NODES=$(NODES) MODE="--block --async" -C parallel $(TEST)
	@cat times.csv

hostfile:
	./make_hostfile.sh $(ROOM) $(SLOTS)
	@cat hostfile

testSave:
	@for v in $(VERSIONS); do $(MAKE) -C $$v testSave; done
	diff ./sequentiel/bin/shalw_512x512_T40.png ./parallel/bin/shalw_512x512_T40.png

clean:
	@for v in $(VERSIONS); do $(MAKE) -C $$v clean; done

mrproper:
	@for v in $(VERSIONS); do $(MAKE) -C $$v mrproper; done
	-rm times.csv hostfile

.PHONY: parallel sequentiel times testSave clean mrproper
