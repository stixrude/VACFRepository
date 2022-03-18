	subroutine flyv(x,n,xm,sig,cor,iflag)

C  Flyvbjerg and Petersen, J. Chem. Phys., 91, 461, 1989 
C  Compute mean value and its uncertainty using proper, non-Gaussian, uncorrelated statistics
C  Added dsig 12/03/14
C
C  Input
C  n = number of data points
C  x = data vector
C  iflag = 0: compute only the mean value, do not compute uncertainty
C  iflag = any other value: compute mean value and uncertainty
C
C  Output
C  xm = mean value of x
C  sig = uncertainty in the mean value of x
C  dsig = uncertainty in sig
C

c	real x(1000000),xp(1000000)
	real x(n),xp(n)
	fn = float(n)

	xm = 0.
	x2m = 0.
	sig = 0.
	do 11 i=1,n
	 xm = xm + x(i)
	 x2m = x2m + x(i)*x(i)
11	continue
	xm = xm/float(n)
	x2m = x2m/float(n)
	if (n .eq. 0) then
	 xm = x(1)
	 x2m = x(1)*x(1)
	end if
	if (iflag .eq. 0 .or. n .le. 2) return

	nrbin = log(fn)/log(2.) - 1
c	print*, 'Number of re-binnings = ',nrbin,xm,n
	ci = (x2m - xm*xm)/float(n-1)
	faci = 1./sqrt(2.*float(n-1))
	dc = faci*ci
	irb = 0
	np = n
c	if (n .lt. 6) write(11,*) 'in flyv',n,xm,x2m,ci,dc,(x(i),i=1,n)

	cmax = -1.
	cold = 0.
c	print 100, irb,np,ci,sqrt(ci),dc
	cmax = ci
	do 3 i=1,n
3	xp(i) = x(i)
	np = fn
	do 1 irb=1,nrbin
	 c = 0.
	 np = np/2
	 do 2 i=1,np
	  if (np .eq. 1) go to 2
	  fac = 1./sqrt(2.*float(np-1))
	  xp(i) = (xp(2*i) + xp(2*i-1))/2
	  c = c + (xp(i) - xm)*(xp(i) - xm)/(float(np)*float(np-1))
c	  if (np .le. 4) print*, irb,np,xp(i),c
2	 continue
	 dc = fac*sqrt(c)
	 diff = sqrt(c) - sqrt(cold)
	 if (abs(diff) .lt. dc) then
c	  print 100, irb,np,c,sqrt(c),dc,'*'
	  if (c .gt. cmax) then
	   cmax = c
	   dcmax = dc
	  end if
	 else
c	  print 100, irb,np,c,sqrt(c),dc
	 end if
	 cold = c
1	continue
	sig = sqrt(cmax)
	dsig = dc
c	print 200, sqrt(cmax),dcmax
c	if (n .lt. 6) write(11,*) 'in flyv',cmax,sig

	cor = float(n)*sig**2/(ci*float(n-1))

	return
100	format(2i10,3f15.6,a5)
200	format(3f15.6,a5)
	end
