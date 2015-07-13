LIBTOOL_DEPS = @LIBTOOL_DEPS@
	libtool: $(LIBTOOL_DEPS)
		$(SHELL) ./config.status --recheck

lib_LTLIBRARIES = libfilter.la
libfilter_la_includedir = $(includedir)/libfilter
libfilter_la_SOURCES = filter.c 
libfilter_la_include_HEADERS = filter.h
libfilter_la_LDFLAGS= -version-info 0:0:0 -release 0.0.1

