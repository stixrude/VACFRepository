	subroutine bfindm(Bg,Ag,fg)
	include 'P1'
	include 'const.inc'
	logical check
        real x(1)
        real fmom0(5)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom0
	ibtyp = 0

C  initial guess 

	if (Bg .lt. 0.) Bg = 0.01
	call bcalc(Ag,Bg,fg)

	x(1) = Bg

c        print '(a20,99e12.5)', 'bfindm initial guess',Ag,Bg,fg

        call broydn(x,1,check)

        Bg = x(1)
        call bcalc(Ag,Bg,fg)

	return
	end
