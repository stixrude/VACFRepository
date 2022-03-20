	function formzfunc(atom)
	include 'valence.inc'

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
