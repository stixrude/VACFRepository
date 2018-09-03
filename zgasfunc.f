	function zgasfunc(f,A,B)
C  Compute gas VDOS Eq. 12 French et al. (2016)
C  Re-arranged following Robert et al. (2015) PRE 91, 033310 Supplementary Material
	logical flag
        double precision rz,zz,u,v
	complex fmempls,fmemmin,imag
	imag = (0.,1.)
	pi = 4.0*atan(1.0)

	zz = 0.
        rz = -pi*f/(sqrt(B))
        call wofz(rz,zz,u,v,flag)
	u = A/2.*sqrt(pi/B)*u
	v = A/2.*sqrt(pi/B)*v
	zgasfunc = 4.*u/(u**2 + v**2 + 4.*pi**2*f**2 + 4.*pi*f*v)

	return
	end
