	subroutine zhunt(z,nfreq,zfind,jlo)
	include 'P1'

	integer jlo,nfreq
	real z(nstepsp),zfind

	jlo = 1
	do 1 i=1,nfreq
	 if (z(i) .lt. zfind) then
	  jlo = i
c	  print*, 'in zhunt',i,z(i),zfind,jlo
	  return
	 end if
1	continue

	return
	end
