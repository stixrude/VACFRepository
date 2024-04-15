	subroutine readxdatcar(natom,ntyp,nstep,atom,wmass,cellmass,density,densitynum,volume,natyp,ityp,acell,tarr,xt)
	include 'P1'
	include 'const.inc'

	character*80 lab,header
        character*80 subs(100)
        character*2 atom(ntypmxp)
	integer natyp(ntypmxp),ityp(natomsp),nchar(100)
	real acell(nstepsp,3,3),wmass(ntypmxp),xt(nstepsp,natomsp,3),tarr(nstepsp),acellalt(nstepsp,3,3)

	open(1,file='XDATCAR',status='old')

	read(1,*) header
        print*, 'header = ',header
        write(12,*) 'header = ',header
        read(1,*) scale
        read(1,*) (acell(1,i,1),i=1,3)
        read(1,*) (acell(1,i,2),i=1,3)
        read(1,*) (acell(1,i,3),i=1,3)
        do 1411 i=1,3
         do 1411 j=1,3
          acell(1,i,j) = scale*acell(1,i,j)
1411    continue
        volume = vcell(acell(1,:,:))
	print*, 'volume = ',volume
	write(12,*) 'volume = ',volume

        read(1,'(a80)') lab
        call parse(lab,subs,nchar,n,80)
        do 142 iatom=1,n
         atom(iatom) = subs(iatom)
         do 1421 iel=1,nelem
          if (atom(iatom) .eq. elnam(iel)) then
           wmass(iatom) = wat(iel)
          end if
1421     continue
142     continue
        print '(1x,a16,i3)', 'number of types ',n
        write(12,'(1x,a16,i3)') 'number of types ',n
        print '(1x,a11,100a3)', 'atom types ',(atom(iatom),iatom=1,n)
        write(12,'(1x,a11,100a3)') 'atom types ',(atom(iatom),iatom=1,n)
        print*, 'atom masses ',(wmass(iatom),iatom=1,n)
        write(12,*) 'atom masses ',(wmass(iatom),iatom=1,n)

        do 112 ntyp=1,ntypmxp
         read(1,*,err=113) (natyp(j),j=1,ntyp)
         backspace 1
112     continue
113     continue
        backspace 1
        ntyp = ntyp - 1
        print*, ntyp,(natyp(j),j=1,ntyp)
        write(12,*) 'atom numbers',(natyp(j),j=1,ntyp)
        wtot = 0.
        do 1422 i=1,ntyp
         wtot = wtot + natyp(i)*wmass(i)
1422    continue
        density = wtot/volume/avn*1.e24
        print*, 'Density = ',density,' g/cm^3'
        write(12,*) 'Density = ',density,' g/cm^3'
        if (n .ne. ntyp) print*, 'inconsistency in atom labels and atom numbers',n,ntyp

        natom = 0
        natpmn = 100000
        do 13 i=1,ntyp
         natom = natom + natyp(i)
         natpmn = min(natyp(i),natpmn)
13      continue
        densitynum = float(natom)/volume
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
21      continue
        print*, 'Cell mass = ',cellmass

        nstep = 0
        do 151 istep=1,nstepsp
         do 143 i=1,3
          do 143 j=1,3
           acell(istep,i,j) = scale*acell(1,i,j)
143      continue
         read(1,*,end=1213) lab
         if (lab .eq. header) then
          read(1,*) scale
          read(1,*) (acell(istep,i,1),i=1,3)
          read(1,*) (acell(istep,i,2),i=1,3)
          read(1,*) (acell(istep,i,3),i=1,3)
ctest          read(1,*) (acellalt(istep,i,1),i=1,3)
ctest          read(1,*) (acellalt(istep,i,2),i=1,3)
ctest          read(1,*) (acellalt(istep,i,3),i=1,3)
          do 141 i=1,3
           do 141 j=1,3
            acell(istep,i,j) = scale*acell(istep,i,j)
ctest           acellalt(istep,i,j) = scale*acellalt(istep,i,j)
141       continue
c          print*, istep,lab,acell(istep,:,:)
          call skip(1,3,ierr)
         end if
         do 121 i=1,natom
          read(1,*,err=1211,end=1213) (xt(istep,i,j),j=1,3)
	  do 122 j=1,3
122	  xt(istep,i,j) = xt(istep,i,j) - floor(xt(istep,i,j))
          go to 1212
1211      print*, 'Error reading from XDATCAR',istep,i,(xt(istep,i,j),j=1,3),(xt(istep,i-1,j),j=1,3)
1212      continue
121      continue
         nstep = nstep + 1
         tarr(istep) = potim*(istep-1)
151     continue
1213    continue
        print*, 'coordinates read in',nstep
        write(12,*) 'coordinates read in',nstep
	do 123 i=1,natom
123	write(129,*) (xt(nstep,i,j),j=1,3)

	return
	end
