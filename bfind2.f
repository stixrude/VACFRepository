	subroutine bfind2(Bg,Ag,fg)
	include 'P1'
	include 'const.inc'
	logical check
	real x(2)
        real fmom0(5)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom0
	ibtyp = 2

C  initial guess 

	call bcalc(Ag,Bg,fg)
	As = Ag
	As = (fmom0(2) - fg*Ag)/(1. - fg)

	x(1) = Bg
	x(2) = As

c        print '(a20,99e12.5)', 'bfind2 initial guess',Ag,Bg,fg,As,fg*Ag+(1.-fg)*As,fmom0(2)

	call broydn(x,1,check)

	Bg = x(1)
	call bcalc(Ag,Bg,fg)

	return
	end
