      SUBROUTINE dftint(func,init,M,NDFT,a,b,w,cosint,sinint)
CLPS  Modify this routine: 
CLPS  Make M and NDFT arguments rather than parameters
CLPS  Assign data(:) an absolute dimension from P1 include file
CLPS  Compute 2pi
CLPS  Make init an argument.  Must be set to a value other than unity for new data set.
	include 'P1'
      INTEGER M,NDFT,MPOL
      REAL a,b,cosint,sinint,w,func,TWOPI
c      PARAMETER (M=64,NDFT=1024,MPOL=6,TWOPI=2.*3.14159265)
c      PARAMETER (MPOL=6,TWOPI=2.*3.14159265)
      PARAMETER (MPOL=6)
      EXTERNAL func
CU    USES dftcor,func,polint,realft
      INTEGER init,j,nn
      REAL aold,bold,c,cdft,cerr,corfac,corim,corre,delta,en,s,sdft,
c     *serr,cpol(MPOL),data(NDFT),endpts(8),spol(MPOL),xpol(MPOL)
     *serr,cpol(MPOL),data(nstepsp),endpts(8),spol(MPOL),xpol(MPOL)
c      SAVE init,aold,bold,delta,data,endpts
c      DATA init/0/,aold/-1.e30/,bold/-1.e30/
      SAVE aold,bold,delta,data,endpts
      DATA aold/-1.e30/,bold/-1.e30/

	TWOPI = 8.0*atan(1.0)
      if (init.ne.1.or.a.ne.aold.or.b.ne.bold.or.data(1).ne.func(a)) then
c        init=1
        aold=a
        bold=b
        delta=(b-a)/M
        do 11 j=1,M+1
          data(j)=func(a+(j-1)*delta)
11      continue
        do 12 j=M+2,NDFT
          data(j)=0.
12      continue
        do 13 j=1,4
          endpts(j)=data(j)
          endpts(j+4)=data(M-3+j)
13      continue
        call realft(data,NDFT,1)
        data(2)=0.
      endif
      en=w*delta*NDFT/TWOPI+1.
      nn=min(max(int(en-0.5*MPOL+1.),1),NDFT/2-MPOL+1)
      do 14 j=1,MPOL
        cpol(j)=data(2*nn-1)
        spol(j)=data(2*nn)
        xpol(j)=nn
        nn=nn+1
14    continue
      call polint(xpol,cpol,MPOL,en,cdft,cerr)
      call polint(xpol,spol,MPOL,en,sdft,serr)
      call dftcor(w,delta,a,b,endpts,corre,corim,corfac)
      cdft=cdft*corfac+corre
      sdft=sdft*corfac+corim
      c=delta*cos(w*a)
      s=delta*sin(w*a)
      cosint=c*cdft-s*sdft
      sinint=s*cdft+c*sdft
      return
      END
