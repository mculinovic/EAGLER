include ../Makefile.rules

MODULE = release

BIN_FILE = $(NAME)
INSTALL_DIR = /usr/local/bin
INSTALLED_BIN = $(INSTALL_DIR)/$(NAME)

CP_FLAGS += -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0


install: all
	@echo [BIN] installing binary to $(INSTALL_DIR)
	@mkdir -p $(INSTALL_DIR)
	@install -s $(BIN_FILE) $(INSTALL_DIR)

uninstall:
	@echo [RM] deleting executable
	@rm -f $(INSTALLED_BIN)
