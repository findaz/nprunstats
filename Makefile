#!/bin/make
ifndef DATE
 $(error DATE must be set!)
endif
UTARGS := $(DATE)

ifndef OUTPUT
 $(error OUTPUT must be set!)
endif
OARGS := $(OUTPUT)

ifdef SXDIR
 SXARGS := --sxdir $(SXDIR)
endif

ifdef CPDIR
 CPARGS := --sxdir $(CPDIR)
endif

ifdef NPROC
 MPARGS := --nproc $(NPROC)
endif

runstats:
	python nprunstats.py $(UTARGS) $(OARGS) $(SXARGS) $(CPARGS) $(MPARGS)

