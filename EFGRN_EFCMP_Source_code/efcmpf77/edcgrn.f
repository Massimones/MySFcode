	subroutine edcgrn(ns,nrec,grndir,grnss0,grnds0,grncl0)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer ns,nrec
	character*80 grndir,grnss0,grnds0,grncl0
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
c	OBSERVATION POSITIONS AND OBSERVABLES
c	=====================================
c
c	(xrec(i),yrec(i),zrec0)=coordinates of the observation positions
c	(Note that zrec0 is fixed)
c	disp = the 3 displcement vector components: ux,uy,uz
c	strain = the 6 strain tensor components: exx,eyy,ezz,exy,eyz,ezx
c	tilt = the two vertical tilt components: dux/dz, duy/dz
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
c	ELASTIC PARAMETERS AT OBSERVATION DEPTH
c	=======================================
c
c	lambda,mu = the two Lame constants in pascal
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision lambda,mu
c
	common/elasticity/lambda,mu
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
c	LOCAL CONSTANTS
c	==============
c	pi, parameter for transforming degree to radian
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision PI,PI2
	data PI,PI2/3.14159265358979d0,6.28318530717959d0/
	double precision DEGTORAD,ZSEPS
	data DEGTORAD,ZSEPS/1.745329252d-02,1.0d-02/
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	GREEN'S FUNNCTIONN PARAMETERS
c	=============================
c
c	Green's function source types:
c	  1 = strike-slip (m12=m21=1)
c	  2 = dip-slip (m13=m31=1)
c	  3 = compensated linear vector dipole (CLVD)
c	      (m11=m22=-1/2, m33=1) (no tangential component)
c	Green's function coordinate system:
c	  (z,r,t) = cylindrical with z being downward(!)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer nr,nz
	double precision r1,r2,z1,z2
	double precision ssdisp(3,NRMAX+2,NZMAX+2)
	double precision ssstrn(6,NRMAX+2,NZMAX+2)
	double precision ssuzr(NRMAX+2,NZMAX+2)
	double precision dsdisp(3,NRMAX+2,NZMAX+2)
	double precision dsstrn(6,NRMAX+2,NZMAX+2)
	double precision dsuzr(NRMAX+2,NZMAX+2)
	double precision cldisp(3,NRMAX+2,NZMAX+2)
	double precision clstrn(6,NRMAX+2,NZMAX+2)
	double precision cluzr(NRMAX+2,NZMAX+2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer i,irec,ips,nps,ir,iz,lendir,nx,ny,idis
	double precision ur,ut,uz,err,ett,ezz,ert,etz,ezr
	double precision x1,x2,y1,y2,zobs,si,co,si2,co2,dis,ddis,azi
	double precision dr,dz,dzs,w00,w10,w01,w11
	double precision uzr,uzt,uzx,uzy,ps,sh,ps2,sh2,oot
	double precision la1,mu1
	character*160 grnss,grnds,grncl
	character*180 dataline
c
c	INITIALIZATION
c	==============
c
	do irec=1,nrec
	  do i=1,3
	    disp(irec,i)=0.d0
	  enddo
	  do i=1,6
	    strain(irec,i)=0.d0
	  enddo
	  do i=1,2
	    tilt(irec,i)=0.d0
	  enddo
	enddo
c
c	READ IN GREEN'S FUNCTIONS
c	=========================
c
	lendir=index(grndir,' ')-1
	grnss=grndir(1:lendir)//grnss0
	grnds=grndir(1:lendir)//grnds0
	grncl=grndir(1:lendir)//grncl0
c
c	READ GREEN'S FUNCTIONS
c	======================
c
c	for point fx source
c
c	ssdisp(1-3): Uz, Ur, Ut
c	ssstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt
c
	print *,'... read Green fcts for the fx source ...'
	open(21,file=grnss,status='old')
	call getdata(21,dataline)
        read(dataline,*)nr,r1,r2,nz,z1,z2,zrec0,lambda,mu
	if(nr.gt.NRMAX)then
	  stop ' Error in edcgrn: to large no of Green fct distantces!'
	else if(nz.gt.NZMAX)then
	  stop ' Error in edcgrn: to large no of Green fct depths!'
	endif
	read(21,*)(((ssdisp(i,ir,iz),i=1,3),(ssstrn(i,ir,iz),i=1,6),
     &            ssuzr(ir,iz),ir=1,nr),iz=1,nz)
	close(21)
c
c	for point fy source
c
c	dsdisp(1-3): Uz, Ur, Ut
c	dsstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt
c
	print *,'... read Green fcts for the fy source ...'
	open(22,file=grnds,status='old')
	call getdata(22,dataline)
        read(dataline,*)nx,x1,x2,ny,y1,y2,zobs,la1,mu1
	if(nx.ne.nr.or.x1.ne.r1.or.x2.ne.r2.or.
     &     ny.ne.nz.or.y1.ne.z1.or.y2.ne.z2)then
	  stop ' Error in edcgrn: different grids in Green functions!'
	endif
	if(la1.ne.lambda.or.mu1.ne.mu)then
	  stop ' Error in edcgrn: different Lame const. in Green fcts!'
	endif
	if(zobs.ne.zrec0)then
	  stop ' Error in edcgrn: different obs.depth in Green fcts!'
	endif
	read(22,*)(((dsdisp(i,ir,iz),i=1,3),(dsstrn(i,ir,iz),i=1,6),
     &            dsuzr(ir,iz),ir=1,nr),iz=1,nz)
	close(22)
c
c	for point fz source
c
c
c	cldisp(1-2): Uz, Ur (Ut=0)
c	clstrn(1-4): Ezz,Err,Ett,Ezr=Erz (Ert=Etr=Etz=Ezt=0)
c
	print *,'... read Green fcts for the fz source ...'
	open(23,file=grncl,status='old')
	call getdata(23,dataline)
        read(dataline,*)nx,x1,x2,ny,y1,y2,zobs,la1,mu1
	if(nx.ne.nr.or.x1.ne.r1.or.x2.ne.r2.or.
     &     ny.ne.nz.or.y1.ne.z1.or.y2.ne.z2)then
	  stop ' Error in edcgrn: different grids in Green functions!'
	endif
	if(la1.ne.lambda.or.mu1.ne.mu)then
	  stop ' Error in edcgrn: different Lame const. in Green fcts!'
	endif
	if(zobs.ne.zrec0)then
	  stop ' Error in edcgrn: different obs.depth in Green fcts!'
	endif
	read(23,*)(((cldisp(i,ir,iz),i=1,3),(clstrn(i,ir,iz),i=1,6),
     &            cluzr(ir,iz),ir=1,nr),iz=1,nz)
	close(23)
c
	do iz=nz+1,nz+2
	  do ir=nr+1,nr+2
	    do i=1,3
	      ssdisp(i,ir,iz)=0.d0
	      dsdisp(i,ir,iz)=0.d0
	      cldisp(i,ir,iz)=0.d0
	    enddo
	    do i=1,6
	      ssstrn(i,ir,iz)=0.d0
	      dsstrn(i,ir,iz)=0.d0
          clstrn(i,ir,iz)=0.d0
	    enddo
	    ssuzr(ir,iz)=0.d0
	    dsuzr(ir,iz)=0.d0
	    cluzr(ir,iz)=0.d0
	  enddo
	enddo
c
	print *,'... read in Green fcts successful ...'
c
c	DISCRETISATION OF RECTANGULAR PLANE SOURCES
c	===========================================
c
c
	dr=(r2-r1)/dble(nr-1)
	dz=(z2-z1)/dble(nz-1)
	if(zrec0.ge.z1.and.zrec0.le.z2.and.r1.eq.0.d0)then
	  iz=idnint((zrec0-z1)/dz)+1
	  dzs=(zrec0-(z1+dz*dble(iz-1)))/dz
	  if(dabs(dzs).le.ZSEPS)then
	    do i=1,3
	      ssdisp(i,1,iz)=ssdisp(i,2,iz)
	      dsdisp(i,1,iz)=dsdisp(i,2,iz)
	      cldisp(i,1,iz)=cldisp(i,2,iz)
	    enddo
	    do i=1,6
	      ssstrn(i,1,iz)=ssstrn(i,2,iz)
	      dsstrn(i,1,iz)=dsstrn(i,2,iz)
	      clstrn(i,1,iz)=clstrn(i,2,iz)
	    enddo
	    ssuzr(1,iz)=ssuzr(2,iz)
	    dsuzr(1,iz)=dsuzr(2,iz)
	    cluzr(1,iz)=cluzr(2,iz)
	  endif
	endif
c
	call edcdisc(ns,nz,z1,z2,dr,dz,nps)
c
c	SUPERPOSITION OF ALL DISCRETE POINT SOURCES
c	===========================================
c
c	disp(1-3): Ux,Uy,Uz
c	strain(1-6): Exx,Eyy,Ezz,Exy=Eyx,Eyz=Ezy,Ezx=Exz
c
	print *,'... superposition of all discrete point sources ...'
	do irec=1,nrec
	  uzx=0.d0
	  uzy=0.d0
	  do ips=1,nps
	    dis=dsqrt((xrec(irec)-pxs(ips))**2
     &               +(yrec(irec)-pys(ips))**2)
	    if(dis.le.0.d0)then
	      azi=0.d0
	    else
	      azi=datan2(yrec(irec)-pys(ips),
     &                   xrec(irec)-pxs(ips))
	    endif
c
	    if(dis.gt.r2)then
	      print *,' Warning: too large distances ignored!'
	      nwarn=nwarn+1
	    else if(pzs(ips).lt.z1)then
	      print *,' Warning: too shallow sources ignored!'
	      nwarn=nwarn+1
	    else if(pzs(ips).gt.z2)then
	      print *,' Warning: too deep sources ignored!'
	      nwarn=nwarn+1
	    else
	      iz=idint((pzs(ips)-z1)/dz)+1
	      dzs=(pzs(ips)-(z1+dz*dble(iz-1)))/dz
	      if(dis.le.r1)then
	        idis=1
	        ddis=0.d0
	      else
	        idis=idint((dis-r1)/dr)+1
	        ddis=(dis-(r1+dr*dble(idis-1)))/dr
	      endif
c
c	      weighting factors for the interpolation
c
	      w00=(1.d0-ddis)*(1.d0-dzs)
	      w10=ddis*(1.d0-dzs)
	      w01=(1.d0-ddis)*dzs
	      w11=ddis*dzs
c
	      co=dcos(azi)
	      si=dsin(azi)
	      co2=dcos(2.d0*azi)
	      si2=dsin(2.d0*azi)
c
c	      pmoment(1-5):
c	        1 = weight for strike-slip: m12=m21=1;
c	        poloidal*sin(2 * theta), toroidal*cos(2 * theta)
c
c	        2 = weight for dip-slip: m13=m31=1
c	        poloidal * cos(theta), toroidal * sin(theta)
c
c	        3 = weight for clvd: m33=-m11=-m22=1
c	        axisymmetric
c
c	        4 = weight for 45 deg strike-slip: m11=-m22=1
c	        greenfct4(theta) = green1(theta + 45 deg)
c
c	        5 = weight for 45 deg dip-slip: m23=m32=1
c	        greenfct5(theta) = green2(theta - 90 deg)
c
c	      contributions from the fx components
c
	      ps=pmoment(1,ips)*si
	      sh=pmoment(1,ips)*co
          oot=pmoment(1,ips)*co2
	      if(pzs(ips).ge.zrec0)then
	        ps2=pmoment(1,ips)*si
	        sh2=pmoment(1,ips)*co
	      else
	        ps2=-pmoment(1,ips)*si
	        sh2=-pmoment(1,ips)*co
	      endif	      
	      uz=ps*(w00*ssdisp(1,idis,iz)+w10*ssdisp(1,idis+1,iz)
     &           +w01*ssdisp(1,idis,iz+1)+w11*ssdisp(1,idis+1,iz+1))
	      ur=ps*(w00*ssdisp(2,idis,iz)+w10*ssdisp(2,idis+1,iz)
     &           +w01*ssdisp(2,idis,iz+1)+w11*ssdisp(2,idis+1,iz+1))
	      ut=sh*(w00*ssdisp(3,idis,iz)+w10*ssdisp(3,idis+1,iz)
     &           +w01*ssdisp(3,idis,iz+1)+w11*ssdisp(3,idis+1,iz+1))
c
	      ezz=ps*(w00*ssstrn(1,idis,iz)+w10*ssstrn(1,idis+1,iz)
     &            +w01*ssstrn(1,idis,iz+1)+w11*ssstrn(1,idis+1,iz+1))
	      err=ps*(w00*ssstrn(2,idis,iz)+w10*ssstrn(2,idis+1,iz)
     &            +w01*ssstrn(2,idis,iz+1)+w11*ssstrn(2,idis+1,iz+1))
	      ett=ps*(w00*ssstrn(3,idis,iz)+w10*ssstrn(3,idis+1,iz)
     &            +w01*ssstrn(3,idis,iz+1)+w11*ssstrn(3,idis+1,iz+1))
	      ezr=ps*(w00*ssstrn(4,idis,iz)+w10*ssstrn(4,idis+1,iz)
     &            +w01*ssstrn(4,idis,iz+1)+w11*ssstrn(4,idis+1,iz+1))
	      ert=sh*(w00*ssstrn(5,idis,iz)+w10*ssstrn(5,idis+1,iz)
     &            +w01*ssstrn(5,idis,iz+1)+w11*ssstrn(5,idis+1,iz+1))
	      etz=sh*(w00*ssstrn(6,idis,iz)+w10*ssstrn(6,idis+1,iz)
     &            +w01*ssstrn(6,idis,iz+1)+w11*ssstrn(6,idis+1,iz+1))
c
	      uzr=ps*(w00*ssuzr(idis,iz)+w10*ssuzr(idis+1,iz)
     &           +w01*ssuzr(idis,iz+1)+w11*ssuzr(idis+1,iz+1))
	      uzt=sh*(w00*ssdisp(1,idis,iz)+w10*ssdisp(1,idis+1,iz)
     &           +w01*ssdisp(1,idis,iz+1)+w11*ssdisp(1,idis+1,iz+1))
     &           *2.d0/dis
c	      contributions from the fy components
c
	      if(pzs(ips).ge.zrec0)then
	        ps2=pmoment(2,ips)*si
	        sh2=pmoment(2,ips)*co
	      else
	        ps2=-pmoment(2,ips)*si
	        sh2=-pmoment(2,ips)*co
	      endif	 
	      ps=pmoment(2,ips)*co
	      sh=pmoment(2,ips)*si
	      ps2=abs(ps)
	      sh2=abs(sh)
	      uz=uz+ps*(w00*dsdisp(1,idis,iz)+w10*dsdisp(1,idis+1,iz)
     &           +w01*dsdisp(1,idis,iz+1)+w11*dsdisp(1,idis+1,iz+1))
	      ur=ur+ps*(w00*dsdisp(2,idis,iz)+w10*dsdisp(2,idis+1,iz)
     &           +w01*dsdisp(2,idis,iz+1)+w11*dsdisp(2,idis+1,iz+1))
	      ut=ut+sh*(w00*dsdisp(3,idis,iz)+w10*dsdisp(3,idis+1,iz)
     &           +w01*dsdisp(3,idis,iz+1)+w11*dsdisp(3,idis+1,iz+1))
c
	      ezz=ezz+ps*(w00*dsstrn(1,idis,iz)+w10*dsstrn(1,idis+1,iz)
     &            +w01*dsstrn(1,idis,iz+1)+w11*dsstrn(1,idis+1,iz+1))
	      err=err+ps*(w00*dsstrn(2,idis,iz)+w10*dsstrn(2,idis+1,iz)
     &            +w01*dsstrn(2,idis,iz+1)+w11*dsstrn(2,idis+1,iz+1))
	      ett=ett+ps*(w00*dsstrn(3,idis,iz)+w10*dsstrn(3,idis+1,iz)
     &            +w01*dsstrn(3,idis,iz+1)+w11*dsstrn(3,idis+1,iz+1))
	      ezr=ezr+ps*(w00*dsstrn(4,idis,iz)+w10*dsstrn(4,idis+1,iz)
     &            +w01*dsstrn(4,idis,iz+1)+w11*dsstrn(4,idis+1,iz+1))
	      ert=ert+sh*(w00*dsstrn(5,idis,iz)+w10*dsstrn(5,idis+1,iz)
     &            +w01*dsstrn(5,idis,iz+1)+w11*dsstrn(5,idis+1,iz+1))
	      etz=etz+sh*(w00*dsstrn(6,idis,iz)+w10*dsstrn(6,idis+1,iz)
     &            +w01*dsstrn(6,idis,iz+1)+w11*dsstrn(6,idis+1,iz+1))
c
	      uzr=uzr+ps*(w00*dsuzr(idis,iz)+w10*dsuzr(idis+1,iz)
     &         +w01*dsuzr(idis,iz+1)+w11*dsuzr(idis+1,iz+1))
	      uzt=uzt+sh*(w00*dsdisp(1,idis,iz)+w10*dsdisp(1,idis+1,iz)
     &           +w01*dsdisp(1,idis,iz+1)+w11*dsdisp(1,idis+1,iz+1))
     &           *(-1.d0/dis)
c
c	      contributions from the fz components
c
	      ps=pmoment(3,ips)
	      sh=abs(pmoment(3,ips))
	      uz=uz+ps*(w00*cldisp(1,idis,iz)+w10*cldisp(1,idis+1,iz)
     &           +w01*cldisp(1,idis,iz+1)+w11*cldisp(1,idis+1,iz+1))
	      ur=ur+ps*(w00*cldisp(2,idis,iz)+w10*cldisp(2,idis+1,iz)
     &           +w01*cldisp(2,idis,iz+1)+w11*cldisp(2,idis+1,iz+1))
	      ut=ut+ps*(w00*cldisp(3,idis,iz)+w10*cldisp(3,idis+1,iz)
     &           +w01*cldisp(3,idis,iz+1)+w11*cldisp(3,idis+1,iz+1))
c
	      ezz=ezz+ps*(w00*clstrn(1,idis,iz)+w10*clstrn(1,idis+1,iz)
     &            +w01*clstrn(1,idis,iz+1)+w11*clstrn(1,idis+1,iz+1))
	      err=err+ps*(w00*clstrn(2,idis,iz)+w10*clstrn(2,idis+1,iz)
     &            +w01*clstrn(2,idis,iz+1)+w11*clstrn(2,idis+1,iz+1))
	      ett=ett+ps*(w00*clstrn(3,idis,iz)+w10*clstrn(3,idis+1,iz)
     &            +w01*clstrn(3,idis,iz+1)+w11*clstrn(3,idis+1,iz+1))
	      ezr=ezr+ps*(w00*clstrn(4,idis,iz)+w10*clstrn(4,idis+1,iz)
     &            +w01*clstrn(4,idis,iz+1)+w11*clstrn(4,idis+1,iz+1))
	      ert=ert+ps*(w00*clstrn(5,idis,iz)+w10*clstrn(5,idis+1,iz)
     &            +w01*clstrn(5,idis,iz+1)+w11*clstrn(5,idis+1,iz+1))
	      etz=etz+ps*(w00*clstrn(6,idis,iz)+w10*clstrn(6,idis+1,iz)
     &            +w01*clstrn(6,idis,iz+1)+w11*clstrn(6,idis+1,iz+1))
c
	      uzr=uzr+ps*(w00*cluzr(idis,iz)+w10*cluzr(idis+1,iz)
     &         +w01*cluzr(idis,iz+1)+w11*cluzr(idis+1,iz+1))
	      uzt=uzt+ps*(w00*cldisp(1,idis,iz)+w10*cldisp(1,idis+1,iz)
     &           +w01*cldisp(1,idis,iz+1)+w11*cldisp(1,idis+1,iz+1))
     &           *(-1.d0/dis)
c
c	      transform to cartesian coordinates
c
	      disp(irec,1)=disp(irec,1)+ur*co-ut*si
	      disp(irec,2)=disp(irec,2)+ur*si+ut*co
	      disp(irec,3)=disp(irec,3)+uz
c
	      strain(irec,1)=strain(irec,1)+err*co*co+ett*si*si-ert*si2
	      strain(irec,2)=strain(irec,2)+err*si*si+ett*co*co+ert*si2
	      strain(irec,3)=strain(irec,3)+ezz
	      strain(irec,4)=strain(irec,4)+0.5d0*(err-ett)*si2+ert*co2
	      strain(irec,5)=strain(irec,5)+ezr*si+etz*co
	      strain(irec,6)=strain(irec,6)+ezr*co-etz*si
c
	      uzx=uzx+uzr*co-uzt*si
	      uzy=uzy+uzr*si+uzt*co
	    endif
	  enddo
c
c	  transform of hrizontal tilts to vertical tilts
c
	  tilt(irec,1)=2.d0*strain(irec,6)-uzx
	  tilt(irec,2)=2.d0*strain(irec,5)-uzy
	enddo
c
	return
	end

