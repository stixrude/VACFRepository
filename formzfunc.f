	function formzfunc(atom)
	parameter (ntable=15)
	character*2 atom
	character*5 atomtab(ntable)
	real z(ntable)
	data ncall/0/
	ncall = ncall + 1
	if (ncall .eq. 1) then
	 open(4,file='/Users/stixrude/work/vasp/codes/sq/formz.txt',status='old')
	 do 1 i=1,ntable
	  read(4,'(a5,5x,9f10.6)') atomtab(i),z(i)
	  write(99,*) atomtab(i),z(i)
1	 continue
	 close (4)
	end if

	do 2 i=1,ntable
	 if (atom .eq. atomtab(i)(1:2)) then
	  iatom = i
	  go to 3
	 end if
2	continue
3	continue

	formzfunc = z(iatom)

	return
	end
