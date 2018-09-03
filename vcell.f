	function vcell(a)

	real a(3,3),gij(3,3),av(3),bv(3),cv(3),vc(3)

        do 12 i=1,3
         do 12 j=1,3
          gij(i,j) = 0.0
          if (i .eq. j) gij(i,j) = 1.0
12       continue

        do 3 i=1,3
         av(i) = a(1,i)
         bv(i) = a(2,i)
         cv(i) = a(3,i)
3       continue

	call cross(av,bv,vc)
        vcell = gdot(vc,cv,gij)

	return
	end

