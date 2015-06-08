include ../Makefile.rules

MODULE = release

BIN_FILE = $(NAME)
INSTALL_DIR = /usr/local/bin
INSTALLED_BIN = $(INSTALL_DIR)/$(NAME)

CXX_FLAGS += -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

INSTALL := install
UNINSTALL := rm
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	# add sudo before install and rm on Linux
	INSTALL := sudo $(INSTALL)
	UNINSTALL := sudo $(UNINSTALL)
endif

install: all
	@echo [BIN] installing binary to $(INSTALL_DIR)
	@mkdir -p $(INSTALL_DIR)
	$(INSTALL) -s $(BIN_FILE) $(INSTALL_DIR)

uninstall:
	@echo [RM] deleting executable
	$(UNINSTALL) -f $(INSTALLED_BIN)
