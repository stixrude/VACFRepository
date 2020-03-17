	subroutine bcalc(Ag,Bg,fg)
        real fmom(5)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,f0,fmom
        pi = 4.0*atan(1.0)

	if (ibtyp .eq. 1) then
	 Ag = 8.*sqrt(Bg/pi)*f0/z0
         fg = f0
	else
         Ag = 4.*Bg/(2. + sqrt(pi*(1. + Bg*z0**2/(4.*g0**(4./5.)*d0**(6./5.)))))
         fg = Ag*z0/8.*sqrt(pi/Bg)
	end if

	return
	end
