	subroutine edgbstab(n)
	use modbes
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer n
c
	include 'edgglobal.h'
c
c 	table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x
c	all multiplied by sqrt(x)
c
	double precision bsdx
	common /bessels/ bsdx
c
	double precision pi,pi2
	data pi,pi2/3.14159265358979d0,6.28318530717959d0/
	integer i,j
	double precision x,xsqrt,a,b
	double precision bessj0,bessj1,bessj
	IF(ALLOCATED(bsfct)) DEALLOCATE(bsfct)
	ALLOCATE(bsfct(0:nnbess1,3))
c
	do j=1,3
	  bsfct(0,j)=0.d0
	enddo
	bsdx=pi2/dble(ndbess)
	  if(n.eq.0)then
c	  horizontal-single-force (fy=F0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj1(x)
	    a=xsqrt*bessj0(x)
	    b=xsqrt*bessj(2,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	else if(n.eq.1)then
c	  horizontal-single-force (fx=F0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj1(x)
	    a=xsqrt*bessj0(x)
	    b=xsqrt*bessj(2,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	else if(n.eq.2)then
c	  vertical-single-force (fz=F0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj0(x)
	    bsfct(i,2)=-xsqrt*bessj1(x)
	    bsfct(i,3)=0.d0
	  enddo
	else if(n.eq.3)then
c	  compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj0(x)
	    bsfct(i,2)=-xsqrt*bessj1(x)
	    bsfct(i,3)=0.d0
	  enddo
	else if(n.eq.4)then
c	  dip-slip (m13=m31=M0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj1(x)
	    a=xsqrt*bessj0(x)
	    b=xsqrt*bessj(2,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	else if(n.eq.5)then
	  write(*,*) 'Debug3'
c	  strike-slip (m12=m21=M0)
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj(2,x)
	    a=xsqrt*bessj1(x)
	    b=xsqrt*bessj(3,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	  write(*,*) 'Debug4'
	else
	  stop ' Error in edgbstab: check 0<= n <= 5?'
	endif
c 
	return
	end
