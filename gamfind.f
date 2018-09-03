	subroutine gamfind(delta,gamma)
	parameter (niter=10,tol=1.e-5,xlow=0.0,xupp=1.0)
	external gameq
	common /gamcom/ deltacom
	deltacom = delta

c	x1 = 0.0
c	x2 = tol
	x1 = 0.5
	x2 = 0.5 + tol
c	print*, 'in gamfind',delta,x1,x2,xlow,xupp
	call cage(gameq,x1,x2,xlow,xupp,ires)
c	print*, 'in gamfind after cage',x1,x2,ires
        gamma = zbrent(gameq,x1,x2,tol)
c	print*, 'in gamfind after zeroin',x1,x2,gamma

	return
	end
