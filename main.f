	include 'P1'
	include 'const.inc'
	character*80 lab,header
	character*80 subs(100)
	character*2 atom(ntypmxp)
	integer natyp(ntypmxp),ityp(natomsp),nchar(100)
	real a(3,3),at(3,3,nstepsp),volt(nstepsp)
	real xt(nstepsp,natomsp,3),wmass(ntypmxp)
	real vt(nstepsp,natomsp,3),vacf(nstepsp,ntypmxp),vacf0(ntypmxp),vacfx0(ntypmxp,3),z(nstepsp,ntypmxp)
	real sigvacf(nstepsp,ntypmxp),vacfint(nstepsp,ntypmxp),vacfmean(nstepsp,ntypmxp)
	real vacfx(nstepsp,ntypmxp,3),zx(nstepsp,ntypmxp,3),dvacf(nstepsp,ntypmxp)
	real vcm(3),sigx(ntypmxp),sigvx(ntypmxp),errav(nstepsp)
	real delta(ntypmxp),gamma(ntypmxp),Aa(ntypmxp),Ba(ntypmxp),fgas(ntypmxp),Wa(ntypmxp)
	real zgas(nstepsp,ntypmxp),zsol(nstepsp,ntypmxp),fmom(5,ntypmxp),fmom0(5),tarr(nstepsp),fmomx(5,ntypmxp)
	real formz(natomsp),diff(ntypmxp),sdiff(ntypmxp),omega0(ntypmxp),omega1(ntypmxp),taumem(ntypmxp)
	real vacfmem(ntypmxp),asec(ntypmxp),bsec(ntypmxp)
	external zfunc,zzfunc
	common /momcom/ ibtyp,ffind,zfind,z0,g0,d0,fmom0
	common /zcom/ jtyp,kxyz,nintegrate,freq,dtime,tarr,vacf,vacfx
	parameter (evk = 11605., elem=1.6021766208e-19)
	parameter (femto=1.e-15)
	parameter (gammax=0.545, iseed=5)
c	parameter (mint=6)
	real cof(mint),xut(mint),xold
	logical, parameter :: solid=.false.
c	real, parameter :: potim=2.0
	pi = 4.0*atan(1.0)
	dtime = potim

	open(3,file='VDATCAR',status='old',iostat=ios)
	if (ios .ne. 0) then
	 print*, 'WARNING no VDATCAR file found.  Opening XDATCAR instead'
	 open(3,file='XDATCAR',status='old')
	end if
	open(12,file='vacfout',status='unknown')
	open(121,file='vacf.txt',status='unknown')
	open(13,file='vdos.txt',status='unknown')
	open(17,file='vgas.txt',status='unknown')

	read(3,*) header
	print*, 'header = ',header
	write(12,*) 'header = ',header
	print*, 'solid =',solid,' potim =',potim
	write(12,*) 'solid =',solid,' potim =',potim
	read(3,*) scale
        read(3,*) (at(i,1,1),i=1,3)
        read(3,*) (at(i,2,1),i=1,3)
        read(3,*) (at(i,3,1),i=1,3)
	do 1411 i=1,3
	 do 1411 j=1,3
	  at(i,j,1) = scale*at(i,j,1)
1411	continue
        do 141 i=1,3
         do 141 j=1,3
          a(i,j) = at(i,j,1)
141     continue
        volt(1) = vcell(at(1,1,1))
        print*, 'volume = ',volt(1)
        write(12,*) 'volume = ',volt(1)
	vol = volt(1)
	volm = vol
	apaco = (vol)**(1./3.)
	read(3,'(a80)') lab
	call parse(lab,subs,nchar,n,80)
	do 142 iatom=1,n
	 atom(iatom) = subs(iatom)
	 do 1421 iel=1,nelem
	  if (atom(iatom) .eq. elnam(iel)) then
	   wmass(iatom) = wat(iel)
	  end if
1421	 continue
142	continue
	print*, 'atom types  ',(atom(iatom),iatom=1,n)
	write(12,*) 'atom types  ',(atom(iatom),iatom=1,n)
	print*, 'atom masses ',(wmass(iatom),iatom=1,n)
	write(12,*) 'atom masses ',(wmass(iatom),iatom=1,n)
        do 112 ntyp=1,ntypmxp
         read(3,*,err=113) (natyp(j),j=1,ntyp)
         backspace 3
112     continue
113     continue
	backspace 3
        ntyp = ntyp - 1
        print*, ntyp,(natyp(j),j=1,ntyp)
        write(12,*) 'atom numbers',(natyp(j),j=1,ntyp)
        wtot = 0.
        do 1422 i=1,ntyp
         wtot = wtot + natyp(i)*wmass(i)
1422    continue
        density = wtot/volt(1)/avn*1.e24
        print*, 'Density = ',density,' g/cm^3'
        write(12,*) 'Density = ',density,' g/cm^3'
	if (n .ne. ntyp) print*, 'inconsistency in atom labels and atom numbers',n,ntyp
	natom = 0
	natpmn = 100000
        do 13 i=1,ntyp
         natom = natom + natyp(i)
	 natpmn = min(natyp(i),natpmn)
13      continue
        densitynum = float(natom)/volt(1)
        print*, 'Number Density = ',densitynum,' atoms/A^3'
        write(12,*) 'Number Density = ',densitynum,' atoms/A^3'
        iityp = 1
        nasum = natyp(iityp)
        do 14 i=1,natom
          if (i .le. nasum) then
           ityp(i) = iityp
          else
           iityp = iityp + 1
           nasum = nasum + natyp(iityp)
           ityp(i) = iityp
          end if
14      continue
	print*, 'total number of atoms = ',natom
	cellmass = 0.
	do 21 i=1,ntyp
	 cellmass = cellmass + natyp(i)*wmass(i)
21	continue
        do 22 iatom=1,natom
         formz(iatom) = formzfunc(atom(ityp(iatom)))
22      continue
	print*, 'Cell mass = ',cellmass
	nstep = 0
	do 151 istep=1,nstepsp
	 read(3,'(a80)',end=1213) lab
	 if (lab .eq. header) then
	  print*, istep,lab
	  call skip(3,7,ierr)
	 end if
	 do 121 i=1,natom
	  read(3,*,err=1211,end=1213) (xt(istep,i,j),j=1,3)
	  go to 1212
1211	  print*, 'Error reading from XDATSUM',istep,i,(xt(istep,i,j),j=1,3),(xt(istep,i-1,j),j=1,3)
1212	  continue
121	 continue
	 nstep = nstep + 1
	 tarr(istep) = potim*(istep-1)
151	continue
1213    continue
	print*, 'coordinates read in',nstep
	ibeg = nstep/ratio
	nseg = nstep - ibeg
c  limit range of integration over time intervals because vacf becomes increasingly noisy at large intervals
	nintegrate = min(nintmax,nseg)
C  Analyze nintegrate in terms of the periodic propagation of sound waves (cf. Haile, "Molecular Dynamics Simulation", Wiley, 1992, pg. 86, eq. 2.115)
	sound = 0.001*(vol)**(1./3.)*1.e-10/(potim*nintegrate*femto)
	bulksound = sound**2*density
	print*, 'Sound speed must be smaller than this value to avoid periodic sound waves = ',sound,' km/s'
	print*, 'Adiabatic bulk modulus must be smaller than this value to avoid periodic sound waves = ',bulksound,' GPa'

C  Compute velocities directly from VDATCAR in which velocities are given in lattice coordinate/time step
C  WARNING: Assumes an orthogonal lattice
	do 134 imom=1,5
	 do 134 jtyp=1,ntyp
	  fmomx(imom,jtyp) = 0.
134	continue
	if (ios .eq. 0) then
	 do 1311 istep=1,nstep
	  time = float(istep-1)*potim
	  do 1312 iatom=1,natom
	   do 1313 j=1,3
            call hunt(tarr,nstep,time,jlo)
            klo = min(max(jlo-(mint-1)/2,1),nstep+1-mint)
	    do 184 mstep=1,mint
	     xut(mstep) = xt(klo+mstep-1,iatom,j)
184	    continue
C  Find coefficients of the interpolating polynomial
            call polcof(tarr(klo),xut(1),mint,cof)
	    vti = 0.
	    vtpi = 0.
	    vtppi = 0.
	    vtpppi = 0.
	    vtppppi = 0.
C  Compute time derivatives
	    do 183 k=2,mint
	     vti = vti + (k-1)*cof(k)*time**(k-2)
	     if (k .gt. 2) vtpi = vtpi + (k-1)*(k-2)*cof(k)*time**(k-3)
	     if (k .gt. 3) vtppi = vtppi + (k-1)*(k-2)*(k-3)*cof(k)*time**(k-4)
	     if (k .gt. 4) vtpppi = vtpppi + (k-1)*(k-2)*(k-3)*(k-4)*cof(k)*time**(k-5)
183	    continue
	    vt(istep,iatom,j) = a(j,j)*xt(istep,iatom,j)/potim
C  Compute moments of the vibrational density of states.  cf. Isbister & McQuarrie (1972) J. Chem. Phys., 56, 736, Eq. 4 and Desjarlais (2013) Eq. 21.
	    if (istep .ge. ibeg) then
	     fmomx(2,ityp(iatom)) = fmomx(2,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vti/potim)**2
	     fmomx(3,ityp(iatom)) = fmomx(3,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtpi/potim)**2
	     fmomx(4,ityp(iatom)) = fmomx(4,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtppi/potim)**2
	     fmomx(5,ityp(iatom)) = fmomx(5,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtpppi/potim)**2
	    end if
1313	   continue
1312	  continue
1311	 continue
	 go to 1314
	end if
C  Compute velocities via polynomial interpolation of atomic positions
C  WARNING: Assumes an orthogonal lattice
	do 131 istep=1,nstep
	 time = float(istep-1)*potim
	 do 132 iatom=1,natom
	  do 133 j=1,3
           call hunt(tarr,nstep,time,jlo)
           klo = min(max(jlo-(mint-1)/2,1),nstep+1-mint)
	   xold = xt(klo,iatom,j)
	   lx = 0
C  Unwrap the portion of the time series needed for polynomial interpolation
	   do 182 mstep=1,mint
	    lx = lx - nint(xt(klo+mstep-1,iatom,j) - xold)
	    xold = xt(klo+mstep-1,iatom,j)
	    xut(mstep) = xt(klo+mstep-1,iatom,j) + lx
182	   continue
C  Find coefficients of the interpolating polynomial
           call polcof(tarr(klo),xut(1),mint,cof)
	   vti = 0.
	   vtpi = 0.
	   vtppi = 0.
	   vtpppi = 0.
	   vtppppi = 0.
C  Compute time derivatives
	   do 181 k=2,mint
	    vti = vti + (k-1)*cof(k)*time**(k-2)
	    if (k .gt. 2) vtpi = vtpi + (k-1)*(k-2)*cof(k)*time**(k-3)
	    if (k .gt. 3) vtppi = vtppi + (k-1)*(k-2)*(k-3)*cof(k)*time**(k-4)
	    if (k .gt. 4) vtpppi = vtpppi + (k-1)*(k-2)*(k-3)*(k-4)*cof(k)*time**(k-5)
	    if (k .gt. 5) vtppppi = vtppppi + (k-1)*(k-2)*(k-3)*(k-4)*(k-5)*cof(k)*time**(k-6)
181	   continue
	   vt(istep,iatom,j) = a(j,j)*vti
C  Compute moments of the vibrational density of states.  cf. Isbister & McQuarrie (1972) J. Chem. Phys., 56, 736, Eq. 4 and Desjarlais (2013) Eq. 21.
	   if (istep .ge. ibeg) then
	    fmomx(2,ityp(iatom)) = fmomx(2,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtpi)**2
	    fmomx(3,ityp(iatom)) = fmomx(3,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtppi)**2
	    fmomx(4,ityp(iatom)) = fmomx(4,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtpppi)**2
	    fmomx(5,ityp(iatom)) = fmomx(5,ityp(iatom)) + wmass(ityp(iatom))*(a(j,j)*vtppppi)**2
	   end if
133	  continue
132	 continue
131	continue
1314	continue

C  Compute Temperature

	avke = 0.
	do 139 j=1,3
139	vcm(j) = 0.
	do 136 istep=ibeg,nstep
         avke1 = 0.
	 do 137 iatom=1,natom
	  do 138 j=1,3
           avke1 = avke1 + 0.5*wmass(ityp(iatom))*vt(istep,iatom,j)**2/float(natom-1)
	   avke = avke + 0.5*wmass(ityp(iatom))*vt(istep,iatom,j)**2/float(natom-1)/float(nseg+1)
	   vcm(j) = vcm(j) + wmass(ityp(iatom))*vt(istep,iatom,j)
138	  continue
137	 continue
         if (istep .gt. ibeg .and. avke1 .gt. 10.*avke*float(nseg+1)/float(istep-ibeg+1)) then
          print*, 'Anomalous velocities',istep,avke*float(nseg+1)/float(istep-ibeg+1),avke1
         end if
136	continue
	do 143 j=1,3
	 vcm(j) = vcm(j)/cellmass/float(nseg+1)
143	continue
	avke = avke*1.e7/(Rgas*evk)
	temp = 2./3.*evk*avke
	cmke = 0.5*cellmass*(vcm(1)**2 + vcm(2)**2 + vcm(3)**2)*1.e7/(Rgas*evk)
	tempcm = 2./3.*evk*cmke
	do 144 ia=1,ntyp
	 fmomx(2,ia) = fmomx(2,ia)*1.e7/(3.*Rgas*temp)/float(nseg+1)/float(natyp(ia))
	 fmomx(3,ia) = fmomx(3,ia)*1.e7/(3.*Rgas*temp)/float(nseg+1)/float(natyp(ia))
	 fmomx(4,ia) = fmomx(4,ia)*1.e7/(3.*Rgas*temp)/float(nseg+1)/float(natyp(ia))
	 fmomx(5,ia) = fmomx(5,ia)*1.e7/(3.*Rgas*temp)/float(nseg+1)/float(natyp(ia))
	 omega0(ia) = sqrt(fmomx(2,ia))
144	continue
	print*, 'Average kinetic energy = ',avke,' eV/atom'
	print*, 'Average temperature = ',temp,' K'
	print*, 'Center of mass velocity, kinetic energy, temperature = ',vcm,cmke,tempcm
	print*, 'Einstein frequencies = ',(omega0(ia)/(2.*pi*femto*1.e12),ia=1,ntyp),' THz'
	print*, 'Einstein frequencies = ',(omega0(ia)/(2.*pi*femto*cspeed),ia=1,ntyp),' cm-1'
        write(12,*) 'Average Temperature = ',temp
	write(12,*) 'Einstein frequencies = ',(omega0(ia)/(2.*pi*femto*1.e12),ia=1,ntyp),' THz'
	write(12,*) 'Einstein frequencies = ',(omega0(ia)/(2.*pi*femto*cspeed),ia=1,ntyp),' cm-1'
	print*, 'Moments from trajectories'
        print '(a21,99f12.5)', '2nd moms(cm-1)',(fmomx(2,jtyp)**(1./2.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        print '(a21,99f12.5)', 'Fourth moments',(fmomx(3,jtyp)**(1./4.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        print '(a21,99f12.5)', 'Sixth moments', (fmomx(4,jtyp)**(1./6.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        print '(a21,99f12.5)', 'Eighth moments',(fmomx(5,jtyp)**(1./8.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,*) 'Moments from trajectories'
        write(12,'(a21,99f12.5)') '2nd moms(cm-1)',(fmomx(2,jtyp)**(1./2.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        write(12,'(a21,99f12.5)') 'Fourth moments',(fmomx(3,jtyp)**(1./4.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        write(12,'(a21,99f12.5)') 'Sixth moments', (fmomx(4,jtyp)**(1./6.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        write(12,'(a21,99f12.5)') 'Eighth moments',(fmomx(5,jtyp)**(1./8.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
        pressureig = densitynum*1.e30*boltzk*temp*1.e-9
        print*, 'Ideal gas pressure = ',pressureig,' GPa'
	if (pressureig .gt. bulksound) print*, 'WARNING: Way over the limit for periodic sound waves.  Try reducing nintegrate.'
	if (pressureig .gt. bulksound) write(12,*) 'WARNING: Way over the limit for periodic sound waves.  Try reducing nintegrate.'

C  Compute Velocity Autocorrelation Function

        do 161 i=1,nstepsp
	 do 163 jtyp=1,ntyp
	  vacf(i,jtyp) = 0.
	  do 164 k=1,3
164	  vacfx(i,jtyp,k) = 0.
163	 continue
161     continue
	iint = 0
	nnint = nseg
	corav = 0.
C  Loop over intervals
	do 15 int=iint,nnint,1
	 if (mod(int,100) .eq. 0) print*, int
C  Loop over time origins
	 do 16 istep=ibeg,nstep-int
         do 165 jtyp=1,ntyp
165      vacfint(istep,jtyp) = 0.
	  do 17 iatom=1,natom
	   do 18 j=1,3
	    tarr(int+1) = int*potim
	    add = vt(istep,iatom,j)*vt(istep+int,iatom,j)
	    vacfint(istep,ityp(iatom)) = vacfint(istep,ityp(iatom)) + add/float(natyp(ityp(iatom)))
	    vacf(int+1,ityp(iatom)) = vacf(int+1,ityp(iatom)) + add/float(nseg+1-int)/float(natyp(ityp(iatom)))
	    vacfx(int+1,ityp(iatom),j) = vacfx(int+1,ityp(iatom),j) + 3.*add/float(nseg+1-int)/float(natyp(ityp(iatom)))
18	   continue
17	  continue
16	 continue
         do 171 jtyp=1,ntyp
          call flyv(vacfint(ibeg,jtyp),nseg+1-int,vacfmean(int+1,jtyp),sigvacf(int+1,jtyp),cor,1)
	  corav = corav + cor/float(ntyp)/float(nseg+1)
171	 continue
15	continue

	do 24 jtyp=1,ntyp
	 vacf0(jtyp) = vacf(iint+1,jtyp)
	 do 2411 k=1,3
2411	 vacfx0(jtyp,k) = vacfx(iint+1,jtyp,k)
	 check = wmass(jtyp)/1000./avn/(3.*boltzk*temp)*vacf0(jtyp)*1.e10
	 check = wmass(jtyp)/(3.*Rgas*temp)*vacf0(jtyp)*1.e7
	 print*, 'check normalization of velocity auto-correlation function',jtyp,check,vacf0(jtyp)
24	continue
	do 241 jtyp=1,ntyp
	 do 242 int=iint,nnint,1
	  vacf(int+1,jtyp) = vacf(int+1,jtyp)/vacf0(jtyp)
	  do 2421 k=1,3
2421	  vacfx(int+1,jtyp,k) = vacfx(int+1,jtyp,k)/vacfx0(jtyp,k)
	  vacfmean(int+1,jtyp) = vacfmean(int+1,jtyp)/vacf0(jtyp)
	  sigvacf(int+1,jtyp) = sigvacf(int+1,jtyp)/vacf0(jtyp)
242	 continue
241	continue
        do 19 int=iint,nnint
	 vacfntyp = 0.
	 do 191 jtyp=1,ntyp
	  sum = 0.
	  sum2 = 0.
	  do 192 k=1,3
	   sum = sum + vacfx(int+1,jtyp,k)/3.
	   sum2 = sum2 + vacfx(int+1,jtyp,k)**2/3.
192	  continue
	  sigvx(jtyp) = sqrt(sum2 - sum**2)
	  vacfntyp = vacfntyp + vacf(int+1,jtyp)/float(ntyp)
191	 continue
C  Uncertainty in VACF computed following ideas on pg. 197 of Allen and Tildesley: Eq. 6.30 modified by 
C  consideration of fewer time origins for longer durations (t_run-t argument near bottom of page).
C  Alternative estimates of uncertainty: variance in vacfx (sigvx), and flyv analysis (sigvacf).
         errav(int+1) = sqrt(2.*corav/float(natom)/float(nseg+1-int)/3.)*(1. - abs(vacfntyp))
         write(121,'(99f21.16)') potim*float(int),(vacf(int+1,jtyp),jtyp=1,ntyp),errav(int+1)
         write(126,'(99f21.16)') potim*float(int),(vacfmean(int+1,jtyp),jtyp=1,ntyp),(sigvacf(int+1,jtyp),jtyp=1,ntyp)
     &                          ,(sigvx(jtyp),jtyp=1,ntyp)
         write(15,'(99f21.16)') potim*float(int),(vacfx(int+1,1,k),k=1,3)
19      continue

C  Cosine Transform

C  Print out frequency in wavenumbers
	df = 1./float(nintegrate)/potim
	nfreq = nintegrate/2.
C  Set up discrete Fourier transform integration dftint (Numerical Recipes, Section 13.9).  
C  M=nintegrate, i.e. the number of time intervals.  
C  N>M and must be a power of 2.  Recommendation is N>4M (comment immediately preceding Eq. 13.9.13 pg. 579.
	ndft = log(float(nintegrate))/log(2.) + 3
	ndft = 2**ndft
	if (ndft .gt. nstepsp) ndft = ndft/2
	print*, 'nintegrate,ndft = ',nintegrate,ndft
C  Gaussian damping prescription from http://cacs.usc.edu/education/phys516/VAC.pdf
c	damp = nintegrate*potim/sqrt(3.)
	damp = nintegrate*potim*dampfactor
c	damp = 1.e15
	call srand(iseed)
	do 223 int=iint,nnint
         time = float(int)*potim
	 do 224 jtyp=1,ntyp
	  dvacf(int+1,jtyp) = abs(vacf(int+1,jtyp))*(1. - exp(-(time/damp)**2))
	  vacf(int+1,jtyp) = vacf(int+1,jtyp)*exp(-(time/damp)**2)
C Test against a model gaussian vacf with the same second moment as that of the MD simulation
C	  vacf(int+1,jtyp) = exp(-0.5*omega0(jtyp)**2*time**2) + errav(int+1)*gasdev(iseed)*exp(-(float(int)*potim/damp)**2)
	  do 226 k=1,3
	   vacfx(int+1,jtyp,k) = vacfx(int+1,jtyp,k)*exp(-(time/damp)**2)
226	  continue
224	 continue
	 write(124,*) time,(vacf(int+1,jtyp),jtyp=1,ntyp),(dvacf(int+1,jtyp),jtyp=1,ntyp)
223	continue
	do 2231 jtyp=1,ntyp
	 if (dvacf(nintegrate,jtyp) .gt. errav(nintegrate)) print*, 'WARNING: Damping exceeds uncertainty.  Consider changing nintegrate.'
	 do 2232 kxyz=1,3
	  init = 0
	  call dftint(zzfunc,init,nintegrate,ndft,0,nintegrate*potim,0.,zx(1,jtyp,kxyz),sinint)
2232	 continue
2231	continue
	kxyz = 0
	do 222 jtyp=1,ntyp
	 do 1 j=1,nfreq
	  f = df*float(j-1)
	  freq = f
          if (mod(j,100) .eq. 0) print*, jtyp,j
	  z(j,jtyp) = 0.
	  init = 1
	  if (j .eq. 1) init = 0
	  call dftint(zzfunc,init,nintegrate,ndft,0,nintegrate*potim,2.*pi*f,z(j,jtyp),sinint)
1	 continue
222	continue
	do 225 j=1,nfreq
	 f = df*float(j-1)
C Compare with model gaussian vacf with the same second moment as that of the MD simulation
	 write(13,*) f/femto*1.e-12,(z(j,jtyp),jtyp=1,ntyp)
     &    ,(sqrt(pi/2.)/omega0(jtyp)*exp(-(2.*pi*f)**2/2./omega0(jtyp)**2),jtyp=1,ntyp)
225	continue

	write(12,*) ' '
C  Estimate uncertainty in self-diffusion coefficient by the difference among the three cartesian directions
	write(12,*) 'Self-diffusion coefficient (m^2/s)'
	do 31 jtyp=1,ntyp
	 sum = 0.
	 sum2 = 0.
	 do 311 k=1,3
	  sum = sum + zx(1,jtyp,k)/3.
	  sum2 = sum2 + zx(1,jtyp,k)*zx(1,jtyp,k)/3.
311	 continue
         diff(jtyp) = vacf0(jtyp)*z(1,jtyp)/1.e5/3.
	 sigx(jtyp) = sqrt(sum2 - sum*sum)
         sdiff(jtyp) = vacf0(jtyp)*sigx(jtyp)/1.e5/3.
C  Compute parameters appearing in the exponential memory function theory of VACF Eq. 7.3.33 Hansen & McDonald, "Theory of Simple Liquids", Academic Press, 2nd Edition, 1986.
	 taumem(jtyp) = Rgas*temp/wmass(jtyp)*1.e3/omega0(jtyp)**2/diff(jtyp)*femto
	 omega1(jtyp) = sqrt(omega0(jtyp)**2 - 1./(2.*taumem(jtyp))**2)
	 tratio = wmass(jtyp)/Rgas/temp*1.e-3*diff(jtyp)*omega0(jtyp)/femto
	 oratio = 1./(2.*omega1(jtyp)*taumem(jtyp))
	 discr = 1. - 4.*omega0(jtyp)**2*taumem(jtyp)**2
         write(12,*) atom(jtyp),diff(jtyp),'+-',sdiff(jtyp)
31	continue

	do 33 int=iint,nnint
         time = float(int)*potim
	 do 34 jtyp=1,ntyp
	  discr = 1. - 4.*omega0(jtyp)**2*taumem(jtyp)**2
	  if (discr .le. 0.) then
	   vacfmem(jtyp) = exp(-time/(2.*taumem(jtyp)))
     &                *(cos(omega1(jtyp)*time)
     &                + 1./(2*omega1(jtyp)*taumem(jtyp))*sin(omega1(jtyp)*time))
	  else
	   ap = 1./2./taumem(jtyp)*(1. - sqrt(discr))
	   am = 1./2./taumem(jtyp)*(1. + sqrt(discr))
	   vacfmem(jtyp) = 1./(ap - am)*(ap*exp(-am*time) - am*exp(-ap*time))
	  end if
34	 continue
         write(1121,'(99f21.16)') time,(vacfmem(jtyp),jtyp=1,ntyp),ap,am,discr
         write(1122,'(99f21.16)') time,(exp(-omega0(jtyp)**2*time**2/2.),jtyp=1,ntyp)
33	continue

        condne = 0.
        do 32 iatom=1,natom
         condne = condne + diff(ityp(iatom))*formz(ityp(iatom))**2
32      continue
        condne = condne*elem*elem/(vol*1.e-30)/(1.3806e-23*temp)
        write(12,*) 'Nernst-Einstein Conductivity = ',condne,' S/m'

C  Calculcate Entropy.  cf. French and Desjarlais (2016) PRE.

C  Normalize z_jtyp.  Eq. (1).
	do 41 i=1,nfreq
	 do 42 jtyp=1,ntyp
	  z(i,jtyp) = 4.*z(i,jtyp)
42	 continue
         f = df*float(i-1)
c	 write(16,*) f/femto/cspeed,(z(i,jtyp),jtyp=1,ntyp)
	 write(16,*) f/femto*1.e-12,(z(i,jtyp)/4.,jtyp=1,ntyp)
41	continue

C  Check sum rule and compute moments.  Eq. 2 and definition of moments in text following Eq. 21 in Desjarlais (2013): M_2n = <omega^(2n)>.
	do 43 jtyp=1,ntyp
	 zfac = 1.0
	 do 431 imom=1,5
431	 fmom(imom,jtyp) = 0.
	 do 43 i=1,nfreq
	  if (z(i,jtyp) .lt. 0.) zfac = 0.0
	  fac = 1.0
	  if (i .eq. 1 .or. i .eq. nfreq) fac = 0.5
	  f = df*float(i-1)
	  do 432 imom=1,5
	   fmom(imom,jtyp) = fmom(imom,jtyp) + zfac*fac*z(i,jtyp)*(2.*pi*f)**(2.*(imom-1))*df
432	  continue
	  write(123,*) f/femto/cspeed,(fmom(imom,1),imom=1,5)
43	continue
	print*, 'Moments from Fourier transform of VACF'
	print*, 'Sum rule check',(fmom(1,jtyp),jtyp=1,ntyp),((fmom(m,jtyp),m=1,5),jtyp=1,ntyp)
	print '(a21,99f12.5)', '2nd moms(cm-1)',(sign(1.,fmom(2,jtyp))*abs(fmom(2,jtyp))**(1./2.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	print '(a21,99f12.5)', 'Fourth moments',(sign(1.,fmom(3,jtyp))*abs(fmom(3,jtyp))**(1./4.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	print '(a21,99f12.5)', 'Sixth moments', (sign(1.,fmom(4,jtyp))*abs(fmom(4,jtyp))**(1./6.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	print '(a21,99f12.5)', 'Eighth moments',(sign(1.,fmom(5,jtyp))*abs(fmom(5,jtyp))**(1./8.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,*) 'Moments from Fourier transform of VACF'
	write(12,'(a21,99f12.5)') '2nd moms(cm-1)',(sign(1.,fmom(2,jtyp))*abs(fmom(2,jtyp))**(1./2.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,'(a21,99f12.5)') 'Fourth moments',(sign(1.,fmom(3,jtyp))*abs(fmom(3,jtyp))**(1./4.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,'(a21,99f12.5)') 'Sixth moments', (sign(1.,fmom(4,jtyp))*abs(fmom(4,jtyp))**(1./6.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,'(a21,99f12.5)') 'Eighth moments',(sign(1.,fmom(5,jtyp))*abs(fmom(5,jtyp))**(1./8.)/(2.*pi*femto*cspeed),jtyp=1,ntyp)

C  Calculate normalized diffusivity, Delta.  Eq. 11.  Assume that N_alpha/V_alpha = N/V, i.e. the one-fluid approximation of Lai et al. (2012) PCCP Eq. 19.
	do 44 jtyp=1,ntyp
	 delta(jtyp) = 2./3.*z(1,jtyp)*femto*sqrt(pi*boltzk*temp/(wmass(jtyp)/1000./avn))*(natom/vol)**(1./3.)*(6./pi)**(2./3.)*1.e10
C  Use moments computed from trajectory
	 if (fmomx(2,jtyp) .gt. 0.) fmom(2,jtyp) = fmomx(2,jtyp)
	 if (fmomx(3,jtyp) .gt. 0.) fmom(3,jtyp) = fmomx(3,jtyp)
	 if (fmomx(4,jtyp) .gt. 0.) fmom(4,jtyp) = fmomx(4,jtyp)
	 if (fmomx(5,jtyp) .gt. 0.) fmom(5,jtyp) = fmomx(5,jtyp)
C  Compute theoretical VACF of the form sech(at)cos(bt) after Isbister and McQuarrie (1972) J. Chem. Phys. 56, 736.
	 C = fmom(3,jtyp)/fmom(2,jtyp)**2
	 asec(jtyp) = 0.5*sqrt((C - 1.)*fmom(2,jtyp))
	 bsec(jtyp) = asec(jtyp)*sqrt((5. - C)/(C - 1.))
	 print*, C,asec(jtyp),bsec(jtyp),asec(jtyp)**2+bsec(jtyp)**2,fmom(2,jtyp)
	 print*, 5.*asec(jtyp)**4 + 6.*asec(jtyp)**2*bsec(jtyp)**2 + bsec(jtyp)**4,fmom(3,jtyp)
	 difsec = Rgas*temp/wmass(jtyp)*1.e3*pi/(2.*asec(jtyp))/cosh(pi*bsec(jtyp)/2.*asec(jtyp))*femto
	 print*, 'sech diffusion coefficients',jtyp,difsec
44	continue
	print*, 'Delta:',(delta(jtyp),jtyp=1,ntyp)
	write(12,'(a6,99f12.5)') 'Delta:',(delta(jtyp),jtyp=1,ntyp)

        do 145 int=iint,nnint
         time = float(int)*potim
	 do 146 jtyp=1,ntyp
	  vacfmem(jtyp) = 1./cosh(asec(jtyp)*time)*cos(bsec(jtyp)*time)
146	 continue
	 write(141,*) time,(vacfmem(jtyp),jtyp=1,ntyp)
145	continue

C  Find gamma by solving numerically Eq. 10.
	do 45 jtyp=1,ntyp
	 call gamfind(delta(jtyp),gamma(jtyp))
45	continue
	print*, 'gamma:',(gamma(jtyp),jtyp=1,ntyp)
	write(12,'(a6,99f12.5)') 'gamma:',(gamma(jtyp),jtyp=1,ntyp)

	entgas2 = 0.
	entgas4 = 0.
	entgasm = 0.
C  Find Ag, Bg, and fg from moments.  For two moment solution (bfind2) use Eqs. 19,20,28 of Desjarlais (2013) PRE.  For four moment solution (bfind4) use Eqs. 19,20,30 of same.
C  Having found fg, calculate W (Eqs. 7,8) and therefore the gas contribution to the entropy (Eq. 5).  For the ideal gas portion, assume that N_alpha/V_alpha = N/V, 
C  i.e. the one-fluid approximation of Lai et al. (2012) PCCP Eq. 19.  Include ideal mixing term: Eq. 13b of Lai et al.
	do 46 jtyp=1,ntyp
	 Ba(jtyp) = fmom(2,jtyp)
	 if (solid) go to 46
	 g0 = gamma(jtyp)
	 d0 = delta(jtyp)
	 z0 = z(1,jtyp)
	 call bcalc(Aa(jtyp),Ba(jtyp),fgas(jtyp))
	 debrog = sqrt(hplanck**2/(2.*pi*wmass(jtyp)/1000./avn*boltzk*temp))
	 xfrac = float(natyp(jtyp))/float(natom)
	 do 461 imom=1,5
461	 fmom0(imom) = fmom(imom,jtyp)
	 if (fmom0(2) .gt. 0. .and. fmom0(3) .gt. 0.) then
	  call bfind2(Ba(jtyp),Aa(jtyp),fgas(jtyp))
	  print '(a15,i5,99e12.5)', 'bfind2',jtyp,Aa(jtyp),Ba(jtyp),fgas(jtyp)
	  Wa(jtyp) = 2.5 - log(debrog**3/(vol*1.e-30)*natyp(jtyp)*fgas(jtyp)) 
     &             + log((1. + gamma(jtyp) + gamma(jtyp)**2 - gamma(jtyp)**3)/(1. - gamma(jtyp))**3) 
     &             + (3.*gamma(jtyp)**2 - 4.*gamma(jtyp))/(1. - gamma(jtyp))**2
	  entgas2 = entgas2 + Wa(jtyp)*fgas(jtyp)*natyp(jtyp)/natom - xfrac*log(xfrac)
	 end if
	 if (fmom0(2) .gt. 0. .and. fmom0(3) .gt. 0. .and. fmom0(4) .gt. 0. .and. fmom0(5) .gt. 0.) then
	  call bfind4(Ba(jtyp),Aa(jtyp),fgas(jtyp))
	  print '(a15,i5,99e12.5)', 'bfind4',jtyp,Aa(jtyp),Ba(jtyp),fgas(jtyp)
	  Wa(jtyp) = 2.5 - log(debrog**3/(vol*1.e-30)*natyp(jtyp)*fgas(jtyp)) 
     &             + log((1. + gamma(jtyp) + gamma(jtyp)**2 - gamma(jtyp)**3)/(1. - gamma(jtyp))**3) 
     &             + (3.*gamma(jtyp)**2 - 4.*gamma(jtyp))/(1. - gamma(jtyp))**2
	  entgas4 = entgas4 + Wa(jtyp)*fgas(jtyp)*natyp(jtyp)/natom - xfrac*log(xfrac)
	 end if
46	continue
C  Find Ag, Bg, and fg so that the high frequency tail of the gas-like component matches the total density of states.  
C  Do this by finding Bg such that f_gas*z_gas(f_match) = z(f_match) where f_match is given by: z(f_match) = z(0)/1000.
C  This follows the argument in Desjarlais (2013) (page 5, pp beginning "Having demonstrated...")
C  The following line skips this step and uses instead the moment matching result to compute zgas.
c	go to 52
	do 51 jtyp=1,ntyp
	 Ba(jtyp) = fmom(2,jtyp)
	 if (solid) go to 51
	 g0 = gamma(jtyp)
	 d0 = delta(jtyp)
	 z0 = z(1,jtyp)
         zfind = z0/1000.
         call hunt(z(1,jtyp),nfreq,zfind,jlo)
	 ffind = df*float(jlo-1)
	 call bfindm(Ba(jtyp),Aa(jtyp),fgas(jtyp))
	 debrog = sqrt(hplanck**2/(2.*pi*wmass(jtyp)/1000./avn*boltzk*temp))
	 xfrac = float(natyp(jtyp))/float(natom)
	 Wa(jtyp) = 2.5 - log(debrog**3/(vol*1.e-30)*natom*fgas(jtyp)) 
     &            + log((1. + gamma(jtyp) + gamma(jtyp)**2 - gamma(jtyp)**3)/(1. - gamma(jtyp))**3) 
     &            + (3.*gamma(jtyp)**2 - 4.*gamma(jtyp))/(1. - gamma(jtyp))**2
	 entgasm = entgasm + Wa(jtyp)*fgas(jtyp)*natyp(jtyp)/natom - xfrac*log(xfrac)
	 print '(a15,i5,99e12.5)', 'bfindm',jtyp,Aa(jtyp),Ba(jtyp),fgas(jtyp),Wa(jtyp),entgasm,gamma(jtyp)
51	continue
52	continue
	print '(a15,99f22.5)', 'sqrt(Ag) (cm-1)',(sqrt(Aa(jtyp))/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,'(a15,99f22.5)') 'sqrt(Ag) (cm-1)',(sqrt(Aa(jtyp))/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	print '(a15,99f22.5)', 'sqrt(Bg) (cm-1)',(sqrt(Ba(jtyp))/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	write(12,'(a15,99f22.5)') 'sqrt(Bg) (cm-1)',(sqrt(Ba(jtyp))/(2.*pi*femto*cspeed),jtyp=1,ntyp)
	print '(a14,99f12.5)', 'gas fractions:',(fgas(jtyp),jtyp=1,ntyp)
	write(12,'(a14,99f12.5)') 'gas fractions:',(fgas(jtyp),jtyp=1,ntyp)
	print '(a3,99f12.5)', 'Wa:',(Wa(jtyp),jtyp=1,ntyp)
	write(12,'(a3,99f12.5)') 'Wa:',(Wa(jtyp),jtyp=1,ntyp)
	print '(a49,3f12.5,a3)', 'Gas entropy (2 moments, 4 moments, tail match) = ',entgas2,entgas4,entgasm,' Nk'
C  Prefer high frequency matching method of computing gas-like VDOS.
	entgas = entgasm

C  Compute gas VDOS Eq. 12 French et al. (2016)
        zz = 0.0
	do 531 imom=1,5
	 do 531 jtyp=1,ntyp
	  fmom(imom,jtyp) = 0.
531	continue
	do 47 i=1,nfreq
         f = df*float(i-1)
	 do 48 jtyp=1,ntyp
	  zgas(i,jtyp) = 0.
	  if (fgas(jtyp) .eq. 0.) go to 48
	  zgas(i,jtyp) = zgasfunc(f,Aa(jtyp),Ba(jtyp))
	  do 532 imom=1,5
	   fmom(imom,jtyp) = fmom(imom,jtyp) + zgas(i,jtyp)*(2.*pi*f)**(2.*(imom-1))*df
532	  continue
48	 continue
c         write(17,*) f/femto/cspeed,(z(i,jtyp),jtyp=1,ntyp),(fgas(jtyp)*zgas(i,jtyp),jtyp=1,ntyp)
         write(17,*) f/femto*1.e-12,(fgas(jtyp)*zgas(i,jtyp)/4.,jtyp=1,ntyp)
47	continue
	print*, 'Test: recalculate moments of gas portion of VDOS'
        print '(a21,99e12.5)', '2nd moms      ',(sign(1.,fmom(2,jtyp))*abs(fmom(2,jtyp))**(1./1.),jtyp=1,ntyp)
     &,Aa(1)
        print '(a21,99e12.5)', 'Fourth moments',(sign(1.,fmom(3,jtyp))*abs(fmom(3,jtyp))**(1./1.),jtyp=1,ntyp)
     &,(Aa(1)**2+2.*Aa(1)*Ba(1))
        print '(a21,99e12.5)', 'Sixth moments', (sign(1.,fmom(4,jtyp))*abs(fmom(4,jtyp))**(1./1.),jtyp=1,ntyp)
     &,(Aa(1)**3 + 4.*Aa(1)**2*Ba(1) + 12.*Aa(1)*Ba(1)**2)
        print '(a21,99e12.5)', 'Eighth moments',(sign(1.,fmom(5,jtyp))*abs(fmom(5,jtyp))**(1./1.),jtyp=1,ntyp)
     &,(Aa(1)**4 + 6.*Aa(1)**3*Ba(1) + 28.*Aa(1)**2*Ba(1)**2 + 120.*Aa(1)*Ba(1)**3)

C  Compute solid entropy Eqs. 4,5,6
	entsol = 0.
	do 49 i=1,nfreq
         f = df*float(i-1)
	 x = hplanck*f/(boltzk*temp)/femto
	 Ws = 0.
	 if (f .gt. 0.) Ws = 3.*(x/(exp(x) - 1.) - log(1. - exp(-x)))
	 write(125,*) f/femto/cspeed,x,Ws,(Wa(jtyp),jtyp=1,ntyp)
	 do 50 jtyp=1,ntyp
	  zsol(i,jtyp) = (z(i,jtyp) - fgas(jtyp)*zgas(i,jtyp))/(1. - fgas(jtyp))
	  entsol = entsol + (1. - fgas(jtyp))*zsol(i,jtyp)*Ws*df*natyp(jtyp)/natom
c	  if (jtyp .eq. 1) write(12,*) f/femto/cspeed,Ws,zsol(i,jtyp),entsol
50	 continue
49	continue
	print*, 'Solid entropy',entsol,' Nk'

C  Compute total entropy
	ent = entgas + entsol
	print*, 'Total entropy',ent,' Nk'
	print*, 'Total entropy',ent/cellmass*Rgas*natom,' J/g/K'
	write (12,'(99a14)') 'Entropy','gas2','gas4','gasm','solid','total(Nk)','total(J/g/K)'
	write(12,'(14x,99f14.5)') entgas2,entgas4,entgasm,entsol,ent,ent/cellmass*Rgas*natom

	stop
	end
