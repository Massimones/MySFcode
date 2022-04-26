	subroutine edcdisc(ns,nz,z1,z2,dr,dz,nps)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer ns,nz,nps
	double precision z1,z2,dr,dz
c
c	inputs:
c	ns = total number of source rectangles
c	nz,z1,z2 = number of depth samples, start and end depths used
c		in Green's functions
c	dlength, dwidth = grid size for discretisation
c
c	returned outputs:
c	nps = total number of discrete point sources
c	other outputs through common blocks
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL CONSTANTS
c	===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision DEGTORAD
	data DEGTORAD/1.745329252E-02/
c
	include 'edcglobal.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	RECTANGULAR SOURCE PLANES
c	=========================
c
c	(xs,ys,zs) = coordinates of the start point of strike
c	with x = north, y = east, z = downward.
c	all angles in degree.
c	NSMAX = the max. number of source rectangles
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision dislocation(NSMAX)
	double precision xs(NSMAX),ys(NSMAX),zs(NSMAX)
	double precision length(NSMAX),width(NSMAX)
	double precision strike(NSMAX),dip(NSMAX),rake(NSMAX)
c
	common/rectangles/dislocation,xs,ys,zs,length,width,
     &                    strike,dip,rake
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	DISTRETE POINT SOURCES
c	======================
c
c	(xs,ys,zs) = coordinates of the discrete point sources
c	with x = north, y = east, z = downward
c	angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision pxs(NPSMAX),pys(NPSMAX),pzs(NPSMAX)
	double precision pmoment(5,NPSMAX)
c
	common/pointsources/pxs,pys,pzs,pmoment
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	OBSERVATION POSITIONS AND OBSERVABLES
c	=====================================
c
c	(xrec(i),yrec(i),zrec0)=coordinates of the observation positions
c	(Note that zrec0 is fixed)
c	disp = the 3 displcement vector components: ux,uy,uz
c	strain = the 6 strain tensor components: exx,eyy,ezz,exy,eyz,ezx
c	NRECMAX = the max. number of observation positions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision xrec(NRECMAX),yrec(NRECMAX)
	double precision zrec0
	double precision disp(NRECMAX,3),strain(NRECMAX,6)
	double precision tilt(NRECMAX,2)
c
	common/obsarray/xrec,yrec,zrec0,disp,strain,tilt
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	WARNING STATISTICS
c	==================
c
c	nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer nwarn
c
	common/warnings/nwarn
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer is,ix,iy,nx,ny,ma1,totlay
	double precision x,y,dx,dy,st,di,ra,magnitude
		double precision dlength,dwidth,smoment
        double precision sm(3,3),disarea,mulay(NPSMAX)
        double precision numlay(200),depthlay(200),mUA(NPSMAX)
        double precision vplay(200),vslay(200),rholay(200)
c cccccccccccc MAX
        smoment=0
c    	open(458,file ='layers_info.txt',status='old')
c    	read(458,*)totlay
c        do ma1=1,totlay
c           read(458,*)numlay(ma1),depthlay(ma1),
c     &           vplay(ma1),vslay(ma1),rholay(ma1)
c        enddo
c          depthlay(totlay)=1000000000
c cccccccccccc MAX
c
	dlength=dr
	nps=0
	do is=1,ns
c
	  st=strike(is)*DEGTORAD
	  di=dip(is)*DEGTORAD
	  ra=rake(is)*DEGTORAD
c
          if(di.gt.0.d0)then
	    dwidth=dmin1(dr,dz/dsin(di))
          else
            dwidth=dr
          endif
c     moment tensor components (double couple)
c	  sm(1,1)=-dsin(di)*dcos(ra)*dsin(2.d0*st)
c     &            -dsin(2.d0*di)*dsin(ra)*(dsin(st))**2
c	  sm(2,2)= dsin(di)*dcos(ra)*sin(2.d0*st)
c     &            -dsin(2.d0*di)*dsin(ra)*(dcos(st))**2
c	  sm(3,3)=-(sm(1,1)+sm(2,2))
c	  sm(1,2)= dsin(di)*dcos(ra)*dcos(2.d0*st)
c     &            +0.5d0*dsin(2.d0*di)*dsin(ra)*dsin(2.d0*st)
c	  sm(2,1)=sm(1,2)
c	  sm(1,3)=-dcos(di)*dcos(ra)*dcos(st)
c     &            -dcos(2.d0*di)*dsin(ra)*dsin(st)
c	  sm(3,1)=sm(1,3)
c	  sm(2,3)=-dcos(di)*dcos(ra)*dsin(st)
c     &            +dcos(2.d0*di)*dsin(ra)*dcos(st)
c	  sm(3,2)=sm(2,3)

c     sm for single forces
	  sm(1,1)=dsin(di)*dcos(st)
	  sm(2,2)= dsin(di)*dsin(st)
	  sm(3,3)=dcos(di)
c
	  nx=max0(1,idnint(length(is)/dlength))
	  ny=max0(1,idnint(width(is)/dwidth))
	  dx=length(is)/dble(nx)
	  dy=width(is)/dble(ny)
c
c	  if one of length and width = 0, then it is a line source
c	  if both length and width = 0, then it is a point source
c
	  disarea=dislocation(is)
	  if(dx.gt.0.d0)then
	    disarea=disarea*dx
	  endif
	  if(dy.gt.0.d0)then
	    disarea=disarea*dy
	  endif
c
	  do ix=1,nx
	    x=dx*(dble(ix)-0.5d0)
	    do iy=1,ny
	      y=dy*(dble(iy)-0.5d0)
	      nps=nps+1
	      if(nps.gt.NPSMAX)then
	        print *,' Warning: too large number for discrete ',
     &                  'point sources (i.e., NPSMAX too small)!'
	        nwarn=nwarn+1
	        nps=NPSMAX
	        return
	      endif
	      pxs(nps)=xs(is)+x*dcos(st)-y*dcos(di)*dsin(st)
	      pys(nps)=ys(is)+x*dsin(st)+y*dcos(di)*dcos(st)
	      pzs(nps)=zs(is)+y*dsin(di)
ccccccccccccc MAX
c          if(dx.gt.0.d0 .AND. dy.gt.0.d0)then
c	        do ma1=1,totlay
c              if(pzs(nps).le.depthlay(ma1)) then
c                mulay(nps)=vslay(ma1)*vslay(ma1)*rholay(ma1)
c                mUA(nps)=mulay(nps)*disarea
c                smoment=smoment+muA(nps)
c                EXIT
c              endif
c            enddo
c               write(459,*)pxs(nps),pys(nps),pzs(nps),mulay(nps),disarea
            
c	      endif
ccccccccccccc MAX
	      if(pzs(nps).lt.z1-dz)then
	        print *,' Warning: parts of source rectangles shallower'
	        print *,'          than the Green function grids!'
	        nwarn=nwarn+1
	      endif
	      if(pzs(nps).gt.z2+dz)then
	        print *,' Warning: parts of source rectangles deeper'
	        print *,'          than the Green function grids!'
	        nwarn=nwarn+1
	      endif
c
c	      1 = weight for strike-slip: m12=m21=1;
c	      2 = weight for dip-slip: m13=m31=1
c	      3 = weight for clvd: m33=-m11=-m22=1
c	      4 = weight for 45 deg strike-slip: m11=-m22=1
c	      5 = weight for 45 deg dip-slip: m23=m32=1
c
c	      pmoment(1,nps)=sm(1,2)*disarea
c	      pmoment(2,nps)=sm(1,3)*disarea
c	      pmoment(3,nps)=sm(3,3)*disarea
c	      pmoment(4,nps)=0.5d0*(sm(1,1)-sm(2,2))*disarea
c	      pmoment(5,nps)=sm(2,3)*disarea
	      pmoment(1,nps)=sm(1,1)*disarea
	      pmoment(2,nps)=sm(2,2)*disarea
	      pmoment(3,nps)=sm(3,3)*disarea
	      pmoment(4,nps)=0.5d0*(sm(1,1)-sm(2,2))*disarea
	      write(*,*) pmoment(1,nps),pmoment(2,nps),pmoment(3,nps)
c
	    enddo
	  enddo
	  write(*,'(a,i6,a,i6,a)')' the ',is,'. components => ',
     &                            nx*ny,' point source'
	enddo
c       magnitude=((2.0/3.0)*log10(smoment*1d7))-10.7
c		write(*,*)'------------------------------------------------'
c		write(*,'(a,i7)')' the total number of point sources: ',nps
c  	 	write(*,*)'---------------M.Nespoli 18/11/16---------------'
c		write(*,*)' the Seismic Moment is: ',smoment, 'Nm'
c   		 write(*,*)' the Mw is: ',magnitude
c
	return
	end

