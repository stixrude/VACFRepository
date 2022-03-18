# Compiles main.f

# CFT     = f77
# LDR     = f77
# mac f77 FFLAGS = -N113 -f -N1 -O
# sgi FFLAGS  = -r8
# CFT     = /Applications/Absoft/bin/f90
# LDR     = /Applications/Absoft/bin/f90
# FFLAGS  = -N113 -lU77 -O -g
CFT = gfortran
LDR = gfortran
FFLAGS  = -ffixed-line-length-136 -fdefault-real-8 -Wall -fno-automatic -O4
COMMAND = main

.f.o :
	$(CFT) $(FFLAGS) $*.f -c

MAIN = main.o

SUBS = \
parse.o skip.o vcell.o cross.o gdot.o flyv.o gamfind.o cage.o zbrent.o gameq.o wofz.o \
broydn.o funcv.o bfind2.o bfind4.o fdjac.o fmin.o lnsrch.o qrdcmp.o qrupdt.o rsolv.o rotate.o \
zfunc.o qromb.o polint.o hunt.o trapzd.o zzfunc.o dftint.o dftcor.o realft.o four1.o \
zgasfunc.o bcalc.o bfindm.o polcof.o formzfunc.o gasdev.o

# LIB1 = /afs/lsa.umich.edu/group/geo/users/stixrude/NUMRECSGI/libsginumrec.a
# LIB1 = /big3/crlb/FORCES/,subroutines/libnumrec.a
# LIB1 = -framework Accelerate
# LIB3 = /Users/stixrude/WORK/SCILIBS/MINPACK/minpack_dp.a
# LIB4 = /Users/stixrude/WORK/SCILIBS/slatec4gf41c/libslatec.a
LIB2 = -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

$(COMMAND): $(MAIN) $(SUBS)
	$(LDR) $(FFLAGS) -o $(COMMAND) $(MAIN) $(SUBS) $(LIB1) $(LIB2) $(LIB3) $(LIB4)
