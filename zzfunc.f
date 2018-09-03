	function zzfunc(t)
	include 'P1'
	real tarr(nstepsp),vacf(nstepsp,ntypmxp)
        real vacfx(nstepsp,ntypmxp,3)
        common /zcom/ jtyp,kxyz,nintegrate,freq,dtime,tarr,vacf,vacfx
	pi = 4.0*atan(1.0)

	itime = t/dtime + 1
	
	if (kxyz .eq. 0) then
	 zzfunc = vacf(itime,jtyp)
	else
	 zzfunc = vacfx(itime,jtyp,kxyz)
	end if

	if (itime .gt. nintegrate+1) zzfunc = 0.

	return
	end
