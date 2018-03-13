VERSIONS:=sequentiel bandes blocks

all:
	@for v in $(VERSIONS); do $(MAKE) -C $$v; done
