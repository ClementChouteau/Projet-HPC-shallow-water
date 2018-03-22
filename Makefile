VERSIONS:=sequentiel parallel

all:
	@for v in $(VERSIONS); do $(MAKE) -C $$v; done
	
clean:
	rm -f obj/*.o

mrproper: clean
	rm -f bin/shalw_par
	rm -f bin/shalw_seq
