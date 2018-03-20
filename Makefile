VERSIONS:=sequentiel parallel

all:
	@for v in $(VERSIONS); do $(MAKE) -C $$v; done