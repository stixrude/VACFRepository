	subroutine bfind4(Bg,Ag,fg)
	include 'P1'
	include 'const.inc'
	logical check
	real x(4)
        real fmom0(5)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom0
	ibtyp = 4

C  initial guess 

        call bcalc(Ag,Bg,fg)
	f1 = (1. - fg)/2.
	A1 = Ag

	x(1) = Bg
	x(2) = A1
	x(3) = f1

c	print '(a20,99e12.5)', 'bfind4 initial guess',Ag,Bg,fg,A1,A2,f1

	call broydn(x,3,check)

	Bg = x(1)
        call bcalc(Ag,Bg,fg)

	return
	end
