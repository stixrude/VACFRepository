	subroutine bcalc(Ag,Bg,fg)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom
        pi = 4.0*atan(1.0)

        Ag = 4.*Bg/(2. + sqrt(pi*(1. + Bg*z0**2/(4.*g0**(4./5.)*d0**(6./5.)))))
        fg = Ag*z0/8.*sqrt(pi/Bg)

	return
	end
