DX = doxygen
DOC = docs/Doxyfile

NAME = eagler

default: release

all: debug release

install:
	@echo [MAKE] install
	@$(MAKE) -C release install

debug:
	@echo [MAKE] $(NAME) $@
	@$(MAKE) -C debug

release:
	@echo [MAKE] $(NAME) $@
	@$(MAKE) -C release

docs:
	@echo [DX] generating documentation
	@$(DX) $(DOC) > /dev/null

clean:
	@echo [MAKE] clean
	@$(MAKE) -C debug clean
	@$(MAKE) -C release clean

uninstall:
	@echo [MAKE] uninstall
	@$(MAKE) -C release uninstall

.PHONY: default all install debug release docs clean uninstall
