	subroutine skip(ifile,nskip,ierr)

	character*1 junk

	ierr = 0
	do 1 i=1,nskip
	 read(ifile,'(a)',end=10,err=10) junk
1	continue
	return

10	continue
	ierr = i-1
	return
	end
