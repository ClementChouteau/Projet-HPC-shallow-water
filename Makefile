VERSIONS:=sequentiel parallel

all:
	@for v in $(VERSIONS); do $(MAKE) -C $$v; done
	
testSave:
	@for v in $(VERSIONS); do $(MAKE) -C $$v testSave; done
	diff ./sequentiel/bin/shalw_512x512_T40.png ./parallel/bin/shalw_512x512_T40.png

clean:
	@for v in $(VERSIONS); do $(MAKE) -C $$v clean; done

mrproper:
	@for v in $(VERSIONS); do $(MAKE) -C $$v mrproper; done
