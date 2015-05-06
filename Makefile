default: all

all: debug release

debug:
	@echo [MAKE] $@
	@$(MAKE) -C debug

release:
	@echo [MAKE] $@
	@$(MAKE) -C release

clean:
	@echo [MAKE] clean
	@$(MAKE) -C debug clean
	@$(MAKE) -C release clean

.PHONY: default all debug release clean
