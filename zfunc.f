	function zfunc(t)
	include 'P1'
	real tarr(nstepsp),vacf(nstepsp,ntypmxp)
        real vacfx(nstepsp,ntypmxp,3)
	integer, parameter :: mintz=7
        common /zcom/ jtyp,kxyz,nintegrate,freq,dtime,tarr,vacf,vacfx
	pi = 4.0*atan(1.0)

	call hunt(tarr,nintegrate+1,t,jlo)
	klo = min(max(jlo-(mintz-1)/2,1),nintegrate+1+1-mintz)
c	print*, t,jlo,klo,tarr(klo),vacf(klo,jtyp),jtyp
	call polint(tarr(klo),vacf(klo,jtyp),mintz,t,vacft,dz)

	zfunc = vacft*cos(2*pi*freq*t)

	return
	end
