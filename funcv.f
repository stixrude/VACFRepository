	subroutine funcv(n,x,fvec)
	real x(n),fvec(n)
	real fmom0(5)
        common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom0
	pi = 4.0*atan(1.0)

	if (ibtyp .eq. 0) then
	 Bg = x(1)
	 call bcalc(Ag,Bg,fg)
         zcalc = fg*zgasfunc(ffind,Ag,Bg)
	 fvec(1) = zcalc - zfind
c	 print*, 'in funcv',ibtyp,Ag,Bg,fg,fvec(1)
	end if

	if (ibtyp .eq. 2) then
	 Bg = x(1)
	 As = x(2)
	 call bcalc(Ag,Bg,fg)
         As = (fmom0(2) - fg*Ag)/(1. - fg)
 	
	 fvec(1) = (1. - fg)*As**2 + fg*(Ag**2 + 2.*Ag*Bg) - fmom0(3)

	 zerotest = (fg*(2.*fmom0(2)*Ag - fmom0(3) - Ag**2) + fmom0(3) - fmom0(2)**2)/(2.*fg*(1. - fg)*Ag) - Bg

c	 print*, 'in funcv',ibtyp,Ag,Bg,fg,As,fmom0(2),fmom0(3),fvec(1),fvec(2),zerotest
	end if

	if (ibtyp .eq. 4) then
	 Bg = x(1)
	 A1 = x(2)
	 f1 = x(3)
	 call bcalc(Ag,Bg,fg)
	 f2 = 1. - f1 - fg
	 A2 = (fmom0(2) - f1*A1 - fg*Ag)/f2

c	 fvec(1) = (f1*A1 +    f2*A2 +    fg*Ag                                                     )**(1./1.) - fmom0(2)**(1./1.)
	 fvec(1) = (f1*A1**2 + f2*A2**2 + fg*(Ag**2 + 2.*Ag*Bg)                                     )**(1./1.) - fmom0(3)**(1./1.)
	 fvec(2) = (f1*A1**3 + f2*A2**3 + fg*(Ag**3 + 4.*Ag**2*Bg + 12.*Ag*Bg**2)                   )**(1./1.) - fmom0(4)**(1./1.)
	 fvec(3) = (f1*A1**4 + f2*A2**4 + fg*(Ag**4 + 6.*Ag**3*Bg + 28.*Ag**2*Bg**2 + 120.*Ag*Bg**3))**(1./1.) - fmom0(5)**(1./1.)

c	 print '(a8,i5,15e13.5)', 'in funcv',ibtyp,Ag,Bg,fg,A1,A2,f1,f2,fmom0(2),fmom0(3),fmom0(4),fmom0(5)
c     &    ,fvec(1),fvec(2),fvec(3),fvec(4)

	end if

	return
	end
