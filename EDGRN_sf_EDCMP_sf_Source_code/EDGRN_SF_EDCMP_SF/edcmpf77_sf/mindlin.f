	subroutine mindlin(ns,NSMAX,nrec,NRECMAX,lambda,mu,dislocations,
     &                 xs,ys,zs,strikes,dips,
     &                 xrec,yrec,zrec0,disp,strain,tilt)
	implicit none
c
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer ns,NSMAX,nrec,NRECMAX
        double precision lambda,mu,nu
        double precision dislocations(NSMAX)
        double precision xs(NSMAX),ys(NSMAX),zs(NSMAX)
        double precision strikes(NSMAX),dips(NSMAX)
        double precision xrec(NRECMAX),yrec(NRECMAX)
        double precision zrec0,ZR,XR,XSO,YSO,YR
        double precision disp(NRECMAX,3),UX,UY,UZ
        double precision strain(NRECMAX,6),tilt(NRECMAX,2)
        double precision pot1,pot2,pot3
        integer j,is,irec
        double precision inclination, AZIMUTH
        double precision EXX,EYY,EZZ,EXY,EXZ,EYZ
        double precision EYX,EZX,EZY
        double precision DEPTH
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       ns = the really used number of rectangular sources
c       NSMAX = the upper limit of ns
c       nrec = the really used number of observation positions
c       NRECMAX = the upper limit of nrec
c
c       lambda, mu = the two Lame constants in Pascal (SI unit)
c
c       (xs,ys,zs) = coordinates of the start point of strike
c       with x = north, y = east, z = downward.
c       all angles in degree.
c       (xrec,yrec,zrec0) = cartesian coordinates of observations
c             Note zrec0 is a fixed constant
c       disp = 3 displacement components: ux,uy,uz
c       strain = 6 strain components: exx,eyy,ezz,exy,eyz,ezx
c       tilt = 2 vertical tilt components: dux/dz, duy/dz
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	from Mindlins's subroutine MIND:
c
	INTEGER IRET
   
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	LOCAL CONSTANTS
c	===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision degtorad,eps
	data degtorad,eps/1.745329252E-02,1.0d-06/ 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision csst,ssst,csdi,ssdi
	double precision disp0(3),tilt0(2),strain0(6)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	PROCESSING
c	==========
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	receiver and source independent variables
c
c
      nu=sngl(lambda/(2*(lambda+mu))) 
	do 901 irec=1,nrec
c	  initialization
c
	  do j=1,6
            strain(irec,j)=0.d0
	  enddo
	  do j=1,3
            disp(irec,j)=0.d0
	  enddo
	  do j=1,2
            tilt(irec,j)=0.d0
	  enddo
c
	  do 900 is=1,ns
c
c
        XR=sngl(xrec(irec))
        YR=sngl(yrec(irec))
        ZR=sngl(zrec0)

	    XSO=sngl(xs(is))
	    YSO=sngl(ys(is))
	    DEPTH=sngl(zs(is))
	    INCLINATION=strikes(is)*degtorad
        AZIMUTH=dips(is)*degtorad
        csst=dcos(INCLINATION)
        ssst=dsin(INCLINATION)
        csdi=dcos(AZIMUTH)
        ssdi=dsin(AZIMUTH)
c
c	      single force
c
	      POT1=sngl(dislocations(is)*csst*ssdi)
	      POT2=sngl(dislocations(is)*ssst*ssdi)
	      POT3=sngl(dislocations(is)*csdi)

	      IRET=1
	      call MIND(XR,YR,ZR,XSO,YSO,DEPTH,INCLINATION,AZIMUTH,
     *             POT1,POT2,POT3,UX,UY,UZ,EXX,EYX,EZX,EXY,
     *             EYY,EZY,EXZ,EYZ,EZZ,IRET,mu,lambda,nu)
c	      if(IRET.eq.1)then
c	        stop ' There is a problem in Okada subroutine!'
c	      endif
c
c	    transform from Okada's to Aki's system
c
            disp0(1)=UX
            disp0(2)=UY
	        disp0(3)=UZ
c
            tilt0(1)=0.D0
            tilt0(2)=0.D0
c
            strain0(1)=EXX
            strain0(2)=EYY
            strain0(3)=EZZ
            strain0(4)=EXY
            strain0(5)=EYZ
            strain0(6)=EZX
c
            do j=1,3
              disp(irec,j)=   disp(irec,j)   + disp0(j)
            enddo
            do j=1,2
              tilt(irec,j)=   tilt(irec,j)   + tilt0(j)
            enddo
            do j=1,6
              strain(irec,j)= strain(irec,j) + strain0(j)
            enddo
900	  continue
901     continue
c
	return
	end
