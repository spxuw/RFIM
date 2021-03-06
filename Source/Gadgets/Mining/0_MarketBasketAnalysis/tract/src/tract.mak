#-----------------------------------------------------------------------
# File    : tract.mak
# Contents: build item and transaction management (on Windows systems)
# Author  : Christian Borgelt
# History : 2008.10.05 file created from apriori makefile
#           2011.05.06 changed to double support reporting/recording
#           2011.08.29 main program fim16 added (mainly for testing)
#           2012.07.27 module tract with write functions added (trawr)
#           2013.04.04 added external modules and tract/train main prgs.
#-----------------------------------------------------------------------
THISDIR  = ..\..\tract\src
UTILDIR  = ..\..\util\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS
CFLAGS   = /nologo /W3 /O2 /GS- $(DEFS) /c $(ADDFLAGS)
INCS     = /I $(UTILDIR)

LD       = link.exe
LDFLAGS  = /nologo /subsystem:console /incremental:no
LIBS     = 

HDRS     = $(UTILDIR)\arrays.h    $(UTILDIR)\memsys.h \
           $(UTILDIR)\symtab.h    $(UTILDIR)\escape.h \
           $(UTILDIR)\tabread.h   $(UTILDIR)\scanner.h \
           tract.h clomax.h report.h fim16.h
OBJS     = $(UTILDIR)\arrays.obj  $(UTILDIR)\memsys.obj \
           $(UTILDIR)\idmap.obj   $(UTILDIR)\escape.obj \
           $(UTILDIR)\tabread.obj $(UTILDIR)\scform.obj \
           tract.obj clomax.obj repcm.obj
PRGS     = fim16 tract train

#-----------------------------------------------------------------------
# Build Module
#-----------------------------------------------------------------------
all:       $(PRGS)

fim16:     $(OBJS) ma16main.obj tract.mak
	$(LD) $(LDFLAGS) $(OBJS) m16main.obj $(LIBS) /out:$@

tract:     $(OBJS) tramain.obj tract.mak
	$(LD) $(LDFLAGS) $(OBJS) tramain.obj $(LIBS) /out:$@

train:     $(OBJS) trnmain.obj $(UTILDIR)\tabwrite.obj tract.mak
	$(LD) $(LDFLAGS) $(OBJS) $(UTILDIR)\tabwrite.obj \
              trnmain.obj $(LIBS) /out:$@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
m16main.obj: $(HDRS)
m16main.obj: fim16.c makefile
	$(CC) $(CFLAGS) $(INCS) fim16.c /Fo$@

tramain.obj: $(HDRS)
tramain.obj: tract.c makefile
	$(CC) $(CFLAGS) $(INCS) tract.c /Fo$@

trnmain.obj: $(HDRS)
trnmain.obj: train.c makefile
	$(CC) $(CFLAGS) $(INCS) train.c /Fo$@

#-----------------------------------------------------------------------
# Frequent Item Set Mining (with at most 16 items)
#-----------------------------------------------------------------------
fim16.obj:   $(HDRS)
fim16.obj:   fim16.c makefile
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN fim16.c /Fo$@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tract.obj:   tract.h $(UTILDIR)\arrays.h  $(UTILDIR)\symtab.h \
                     $(UTILDIR)\tabread.h
tract.obj:   tract.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN tract.c /Fo$@

trawr.obj:   tract.h $(UTILDIR)\arrays.h  $(UTILDIR)\symtab.h \
                     $(UTILDIR)\tabread.h $(UTILDIR)\tabwrite.h
trawr.obj:   tract.c makefile
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN /D TA_WRITE tract.c /Fo$@

tatree.obj: tract.h $(UTILDIR)\arrays.h $(UTILDIR)\symtab.h
tatree.obj: tract.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN /D TATREEFN tract.c /Fo$@

#-----------------------------------------------------------------------
# Train Management
#-----------------------------------------------------------------------
train.obj: train.h $(UTILDIR)\arrays.h  $(UTILDIR)\symtab.h \
                   $(UTILDIR)\tabread.h
train.obj: train.c makefile
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN train.c /Fo$@

trnwr.obj: train.h $(UTILDIR)\arrays.h  $(UTILDIR)\symtab.h \
                   $(UTILDIR)\tabread.h $(UTILDIR)\tabwrite.h
trnwr.obj: train.c makefile
	$(CC) $(CFLAGS) $(INCS) /D NOMAIN /D TRN_WRITE train.c /Fo$@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
report.obj: report.h tract.h $(UTILDIR)\symtab.h
report.obj: report.c tract.mak
	$(CC) $(CFLAGS) $(INCS) report.c /Fo$@

repdbl.obj: report.h tract.h $(UTILDIR)\symtab.h
repdbl.obj: report.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D SUPP_T=double report.c /Fo$@

repcm.obj:  report.h tract.h $(UTILDIR)/symtab.h
repcm.obj:  report.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D ISR_CLOMAX report.c /Fo$@

repcmd.obj: report.h tract.h $(UTILDIR)/symtab.h
repcmd.obj: report.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D SUPP_T=double /D ISR_CLOMAX report.c /Fo$@

#-----------------------------------------------------------------------
# Prefix Tree Management for Closed/Maximal Filtering
#-----------------------------------------------------------------------
clomax.obj: clomax.h $(UTILDIR)\arrays.h $(UTILDIR)\memsys.h
clomax.obj: clomax.c tract.mak
	$(CC) $(CFLAGS) $(INCS) -c clomax.c /Fo$@

cmdbl.obj:  clomax.h $(UTILDIR)\arrays.h $(UTILDIR)\memsys.h
cmdbl.obj:  clomax.c tract.mak
	$(CC) $(CFLAGS) $(INCS) /D SUPP_T=double -c clomax.c /Fo$@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)\arrays.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak arrays.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\memsys.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak memsys.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\idmap.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak idmap.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\escape.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak escape.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\tabread.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak tabread.obj ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\scform.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak scform.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	-@erase /Q *~ *.obj *.idb *.pch $(PRGS)
