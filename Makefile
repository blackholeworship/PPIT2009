########################################################################
# Generic Makefile for RNPL 'f77' application.
#
# Requires following environment variables to be set to appropriate
# values:
#
# RNPL_RNPL 
# RNPL_RNPL_FLAGS
# RNPL_F77 
# RNPL_F77LOAD
# RNPL_F77PP 
# RNPL_FLIBS
#
########################################################################
.IGNORE:

########################################################################
#  NOTE: This Makefile uses the Bourne shell, 'sh'
########################################################################
SHELL = /bin/sh

########################################################################
#  Set 'APP' to application name stem (prefix) then execute
#  'make fix' to convert Makefile to use explicit targets
########################################################################
APP        =  emkg

########################################################################
# If your application uses headers and/or libraries from 
# non-system locations, define the following macros appropriately ...
# (set to white-space separated path names, don't use 'csh' ~ notation
# for home directories)
########################################################################
USER_INC_PATHS =
USER_LIB_PATHS =

########################################################################
# If you want to set non-default flags for the 'f77' compiler, do so
# here
########################################################################
RNPL      = $(RNPL_RNPL)

F77       = $(RNPL_F77)
F77_LOAD  = $(RNPL_F77LOAD)
F77PP     = $(RNPL_F77PP)
FLIBS     = $(RNPL_FLIBS)

.f.o:
	$(F77) -c $*.f

all: $(APP) $(APP)_init

fix: Makefile
	sed "s@$(APP)@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

$(APP).f:      $(APP)_rnpl elliptics_init.inc elliptics.inc
	$(RNPL) -l allf   $(RNPL_RNPL_FLAGS) $(APP)_rnpl

updates.f:       $(APP)_rnpl elliptics_init.inc elliptics.inc
initializers.f:  $(APP)_rnpl elliptics_init.inc elliptics.inc
residuals.f:     $(APP)_rnpl elliptics_init.inc elliptics.inc
$(APP)_init.f: $(APP)_rnpl elliptics_init.inc elliptics.inc

$(APP)_init: $(APP)_init.o updates.o initializers.o residuals.o 
	$(F77_LOAD) $(APP)_init.o updates.o initializers.o $ residuals.o $(FLIBS) -o $(APP)_init

$(APP): $(APP).o updates.o residuals.o elliptics_init.inc elliptics.inc $(APP)_rnpl
	$(F77_LOAD) $(APP).o updates.o residuals.o $(FLIBS) -o $(APP)

clean: 
	 rm *.exe
	 rm make.out do.out > /dev/null 2>&1
	 rm *.hdf *.sdf .rn*  $(APP)  $(APP)_init > /dev/null 2>&1
	 rm *.o > /dev/null 2>&1
	 rm initfrag.f updates.f initializers.f residuals.f $(APP).f $(APP)_init.f 
	 rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc 
	 rm *.dump
