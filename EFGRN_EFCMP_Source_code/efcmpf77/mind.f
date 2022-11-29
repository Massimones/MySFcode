	  SUBROUTINE MIND(XR,YR,ZR,XSO,YSO,DEPTH,INCLINATION,AZIMUTH,
     *             POT1,POT2,POT3,UX,UY,UZ,EXX,EYX,EZX,EXY,
     *             EYY,EZY,EXZ,EYZ,EZZ,IRET,mu,lambda,nu)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      double precision   XR,YR,ZR,DEPTH,DIP,POT1,POT2,POT3,DU,               
     *         UX,UY,UZ,EXX,EYX,EZX,EXY,EYY,EZY,EXZ,EYZ,EZZ,
     *         R1,R2,r,XP,YP,ZP,mu,lambda,nu,
     *         SXX,SYX,SZX,SYY,SZY,SZZ,SR,ST,SRZ,
     *         matrScar,matrScyl,matrRot,INCLINATION
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	LOCAL CONSTANTS
c	===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DIMENSION  U(12),DU(12),matrScar(3,3),matrScyl(3,3),matrRot(3,3)                       
      DATA  F0/0.D0/   
      DATA pi/3.141592653589793D0/                                                     
      UX=F0                                                         
      UY=F0                                                             
      UZ=F0                                                             
      EXX=F0                                                            
      EYX=F0                                                            
      EZX=F0                                                            
      EXY=F0                                                            
      EYY=F0                                                            
      EZY=F0                                                            
      EXZ=F0                                                            
      EYZ=F0                                                            
      EZZ=F0  
      XP=XR-XSO
      YP=YR-YSO
      ZP=ZR-DEPTH
      r= SQRT(XP**2 + YP**2)    
      R1=SQRT(r**2+(ZR-DEPTH)**2)
      R2=SQRT(r**2+(ZR+DEPTH)**2)                                             
C=====  FX CONTRIBUTION  ==============                                 
C======================================
       DO  I=1,12 
        U(I)=F0  
       END DO                      
      IF(POT1.NE.F0) THEN   
        DU( 1)=XP*YP/(16*pi*mu*(1-nu))*(1/(R1**3)
     &    +(3-4*nu)/(R2**3)
     &    -6*DEPTH*ZR/(R2**5)
     &    -(4*(1-nu)*(1-2*nu))/(R2*(R2+DEPTH+ZR)**2))                                            
        DU( 2)=1/(16*pi*mu*(1-nu))*((3-4*nu)/R1+1/R2+(YP**2)/(R1**3)
     &   +(3-4*nu)*(YP**2)/(R2**3)
     &   +(2*DEPTH*ZR)/(R2**3)*(1-(3*(YP**2)/(R2**2)))
     &   +4*(1-nu)*(1-2*nu)/(R2+DEPTH+ZR)*(1
     &   -(YP**2)/(R2*(R2+DEPTH+ZR))))                                            
        DU( 3)= YP/(16*pi*mu*(1-nu))*((ZR-DEPTH)/(R1**3)
     &    +(3-4*nu)*(ZR-DEPTH)/(R2**3)-6*DEPTH*ZR*(ZR+DEPTH)/(R2**5)
     &    +4*(1-nu)*(1-2*nu)/(R2*(R2+ZR+DEPTH)))                                             
        SYY=YP/(8*pi*(1-nu))*(-(1-2*nu)/R1**3
     &    +(1-2*nu)*(5-4*nu)/R2**3
     &    -3*YP**2/R1**5-3*((3-4*nu)*YP**2)/R2**5
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(3
     &    -YP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &    +6*(DEPTH/R2**5)*(3*DEPTH-(3-2*nu)*(DEPTH+ZR)
     &    +5*YP**2*ZR/R2**2))                            
        SYX= XP/(8*pi*(1-nu))*(-(1-2*nu)/R1**3
     &    +(1-2*nu)/R2**3-3*YP**2/R1**5
     &   -3*(3-4*nu)*YP**2/R2**5
     &   -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(1
     &   -YP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &   -6*DEPTH*ZR/R2**5*(1-5*YP**2/R2**2))                      
        SZY=  1/(8*pi*(1-nu))*(-(1-2*nu)*(ZR-DEPTH)/R1**3
     &   +(1-2*nu)*(ZR-DEPTH)/R2**5-3*YP**2*(ZR-DEPTH)/R1**5
     &   -3*(3-4*nu)*YP**2*(DEPTH+ZR)/R2**5
     &   -6*DEPTH/R2**5*(ZR*(ZR+DEPTH)-(1-2*nu)*YP**2
     &   -5*YP**2*ZR*(ZR+DEPTH)/R2**2))                                             
        SXX= YP/(8*pi*(1-nu))*((1-2*nu)/R1**3
     &    +(1-2*nu)*(3-4*nu)/R2**3
     &    -3*XP**2/R1**5-3*((3-4*nu)*XP**2)/R2**5
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(1
     &    -XP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &    +6*(DEPTH/R2**5)*(DEPTH-(1-2*nu)*(DEPTH+ZR)
     &    +5*XP**2*ZR/R2**2))                    
        SZX= YP*XP/(8*pi*(1-nu))*(-3*(ZR-DEPTH)/R1**5
     &    -3*(3-4*nu)*(ZR+DEPTH)/R2**5
     &    +6*DEPTH/R2**5*(1-2*nu+5*ZR*(ZR+DEPTH)/R2**2))                                          
        SZZ= YP/(8*pi*(1-nu))*((1-2*nu)/R1**3
     &    -(1-2*nu)/R2**3
     &    -3*(DEPTH-ZR)**2/R1**5
     &    -3*((3-4*nu)*(DEPTH+ZR)**2)/R2**5
     &    +6*(DEPTH/R2**5)*(DEPTH+(1-2*nu)*(DEPTH+ZR)
     &    +5*ZR*(DEPTH+ZR)**2/R2**2))   
        DU( 4)= 1/(2*mu)*(SXX-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                           
        DU( 5)= 1/(2*mu)*SYX                       
        DU( 6)= 1/(2*mu)*SZX                             
        DU( 7)= 1/(2*mu)*SYX                      
        DU( 8)= 1/(2*mu)*(SYY-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                    
        DU( 9)= 1/(2*mu)*SZY                     
        DU(10)= 1/(2*mu)*SZX                     
        DU(11)= 1/(2*mu)*SZY                        
        DU(12)= 1/(2*mu)*(SZZ-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))  
        DO I=1,12                                                   
        U(I)=U(I)+POT1*DU(I) 
        END DO                                        
      ENDIF      
      IF(POT2.NE.F0) THEN                                               
        DU( 2)=XP*YP/(16*pi*mu*(1-nu))*(1/(R1**3)
     &    +(3-4*nu)/(R2**3)
     &    -6*DEPTH*ZR/(R2**5)
     &    -(4*(1-nu)*(1-2*nu))/(R2*(R2+DEPTH+ZR)**2))                                            
        DU( 1)=1/(16*pi*mu*(1-nu))*((3-4*nu)/R1+1/R2+(XP**2)/(R1**3)
     &   +(3-4*nu)*(XP**2)/(R2**3)
     &   +(2*DEPTH*ZR)/(R2**3)*(1-(3*(XP**2)/(R2**2)))
     &   +4*(1-nu)*(1-2*nu)/(R2+DEPTH+ZR)*(1
     &   -(XP**2)/(R2*(R2+DEPTH+ZR))))                                            
        DU( 3)= XP/(16*pi*mu*(1-nu))*((ZR-DEPTH)/(R1**3)
     &    +(3-4*nu)*(ZR-DEPTH)/(R2**3)-6*DEPTH*ZR*(ZR+DEPTH)/(R2**5)
     &    +4*(1-nu)*(1-2*nu)/(R2*(R2+ZR+DEPTH)))                  
        SXX=XP/(8*pi*(1-nu))*(-(1-2*nu)/R1**3
     &    +(1-2*nu)*(5-4*nu)/R2**3
     &    -3*XP**2/R1**5-3*((3-4*nu)*XP**2)/R2**5
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(3
     &    -XP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &    +6*(DEPTH/R2**5)*(3*DEPTH-(3-2*nu)*(DEPTH+ZR)
     &    +5*XP**2*ZR/R2**2))                            
        SYX= YP/(8*pi*(1-nu))*(-(1-2*nu)/R1**3
     &    +(1-2*nu)/R2**3-3*XP**2/R1**5
     &   -3*(3-4*nu)*XP**2/R2**5
     &   -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(1
     &   -XP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &   -6*DEPTH*ZR/R2**5*(1-5*XP**2/R2**2))                      
        SZX=  1/(8*pi*(1-nu))*(-(1-2*nu)*(ZR-DEPTH)/R1**3
     &   +(1-2*nu)*(ZR-DEPTH)/R2**5-3*XP**2*(ZR-DEPTH)/R1**5
     &   -3*(3-4*nu)*XP**2*(DEPTH+ZR)/R2**5
     &   -6*DEPTH/R2**5*(ZR*(ZR+DEPTH)-(1-2*nu)*XP**2
     &   -5*XP**2*ZR*(ZR+DEPTH)/R2**2))                                             
        SYY= XP/(8*pi*(1-nu))*((1-2*nu)/R1**3
     &    +(1-2*nu)*(3-4*nu)/R2**3
     &    -3*YP**2/R1**5-3*((3-4*nu)*YP**2)/R2**5
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+DEPTH+ZR)**2)*(1
     &    -YP**2*(3*R2+DEPTH+ZR)/(R2**2*(R2+DEPTH+ZR)))
     &    +6*(DEPTH/R2**5)*(DEPTH-(1-2*nu)*(DEPTH+ZR)
     &    +5*YP**2*ZR/R2**2))                    
        SZY= XP*YP/(8*pi*(1-nu))*(-3*(ZR-DEPTH)/R1**5
     &    -3*(3-4*nu)*(ZR+DEPTH)/R2**5
     &    +6*DEPTH/R2**5*(1-2*nu+5*ZR*(ZR+DEPTH)/R2**2))                                          
        SZZ= XP/(8*pi*(1-nu))*((1-2*nu)/R1**3
     &    -(1-2*nu)/R2**3
     &    -3*(DEPTH-ZR)**2/R1**5
     &    -3*((3-4*nu)*(DEPTH+ZR)**2)/R2**5
     &    +6*(DEPTH/R2**5)*(DEPTH+(1-2*nu)*(DEPTH+ZR)
     &    +5*ZR*(DEPTH+ZR)**2/R2**2))   
        DU( 4)= 1/(2*mu)*(SXX-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                           
        DU( 5)= 1/(2*mu)*SYX                       
        DU( 6)= 1/(2*mu)*SZX                             
        DU( 7)= 1/(2*mu)*SYX                      
        DU( 8)= 1/(2*mu)*(SYY-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                    
        DU( 9)= 1/(2*mu)*SZY                     
        DU(10)= 1/(2*mu)*SZX                     
        DU(11)= 1/(2*mu)*SZY                        
        DU(12)= 1/(2*mu)*(SZZ-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))           
        DO I=1,12                                                   
        U(I)=U(I)+POT2*DU(I)   
        END DO                                     
      ENDIF   
      IF(POT3.NE.F0) THEN     
      IF(r.NE.F0) THEN                                        
        DU( 2)= (r/(16*pi*mu*(1-nu))*((ZR-DEPTH)/(R1**3)
     &    +(3-4*nu)*(ZR-DEPTH)/(R2**3)
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+ZR+DEPTH))
     &    +6*DEPTH*ZR*(ZR+DEPTH)/(R2**5)))*(YP/r)  
        DU( 1)= (r/(16*pi*mu*(1-nu))*((ZR-DEPTH)/(R1**3)
     &    +(3-4*nu)*(ZR-DEPTH)/(R2**3)
     &    -4*(1-nu)*(1-2*nu)/(R2*(R2+ZR+DEPTH))
     &    +6*DEPTH*ZR*(ZR+DEPTH)/(R2**5)))*(XP/r)   
      ELSE  
        DU( 1)=F0 
        DU( 2)=F0    
      ENDIF                  
        DU( 3)= 1/(16*pi*mu*(1-nu))*((3-4*nu)/(R1)
     &    +(8*(1-nu)**2-(3-4*nu))/R2
     &    +(ZR-DEPTH)**2/(R1**3)
     &    +((3-4*nu)*(ZR+DEPTH)**2-2*ZR*DEPTH)/(R2**3)
     &    +(6*ZR*DEPTH*(ZR+DEPTH)**2)/(R2**5)) 
        SR=1/(8*pi*(1-nu))*((1-2*nu)*(ZR-DEPTH)/R1**3
     &    -(1-2*nu)*(ZR+7*DEPTH)/R2**3+4*(1-nu)*(1
     &    -2*nu)/(R2*(R2+ZR+DEPTH))-3*r**2*(ZR-DEPTH)/R1**5
     &    +(6*DEPTH*(1-2*nu)*(ZR+DEPTH)**2-6*DEPTH**2*(ZR+DEPTH)
     &    -3*(3-4*nu)*r**2*(ZR-DEPTH))/R2**5
     &    -30*DEPTH*r**2*ZR*(ZR+DEPTH)/R2**7)
        ST=(1-2*nu)/(8*pi*(1-nu))*((ZR-DEPTH)/R1**3
     &    +((3-4*nu)*(ZR+DEPTH)-6*DEPTH)/R2**3
     &    -4*(1-nu)/(R2*(R2+DEPTH+ZR))
     &    +6*DEPTH*(DEPTH+ZR)**2/R2**5
     &    -6*DEPTH**2*(ZR+DEPTH)/((1-2*nu)*R2**5))
        SZZ= 1/(8*pi*(1-nu))*(-(1-2*nu)*(ZR-DEPTH)/R1**3
     &    +(1-2*nu)*(ZR-DEPTH)/R2**3-3*(ZR-DEPTH)**3/R1**5
     &    -(3*(3-4*nu)*ZR*(ZR+DEPTH)**2
     &    -3*DEPTH*(ZR+DEPTH)*(5*ZR-DEPTH))/R2**5
     &    -30*ZR*DEPTH*(ZR+DEPTH)**3/R2**7)
        SRZ=r/(8*pi*(1-nu))*(-(1-2*nu)/R1**3+(1-2*nu)/R2**3
     &    -3*(ZR-DEPTH)**2/R1**5-(3*(3-4*nu)*ZR*(ZR+DEPTH)
     &    -3*DEPTH*(3*ZR+DEPTH))/R2**5
     &    -30*DEPTH*ZR*(ZR+DEPTH)**2/R2**7)
      IF(YP.NE.F0) THEN 
        matrRot(1,1)= YP/r   
        matrRot(2,2)= YP/r
      ELSE
        matrRot(1,1)= 1   
        matrRot(2,2)= 1
      ENDIF      
      IF(XP.NE.F0) THEN 
        matrRot(1,2)= -XP/r   
        matrRot(2,1)= XP/r
      ELSE
        matrRot(1,2)= -1   
        matrRot(2,1)= 1
      ENDIF     
        matrRot(1,3)= 0.d0
        matrRot(2,3)= 0
        matrRot(3,1)= 0
        matrRot(3,2)= 0
        matrRot(3,3)= 1
        matrScyl(1,1)= SR      
        matrScyl(1,2)= 0.d0
        matrScyl(1,3)= SRZ
        matrScyl(2,1)= 0.d0
        matrScyl(2,2)= ST
        matrScyl(2,3)= 0.d0
        matrScyl(3,1)= SRZ
        matrScyl(3,2)= 0.d0
        matrScyl(3,3)= SZZ
        matrScar=matmul(matrRot,matrScyl)
        matrScar=matmul(matrScar,transpose(matrRot))
        SYY= matrScar(1,1)                          
        SYX= matrScar(1,2)                    
        SZY= matrScar(1,3)                                       
        SXX= matrScar(2,2)                 
        SZX= matrScar(2,3)
        SZZ= matrScar(3,3)                                      
        DU( 4)= 1/(2*mu)*(SXX-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                           
        DU( 5)= 1/(2*mu)*SYX                       
        DU( 6)= 1/(2*mu)*SZX                             
        DU( 7)= 1/(2*mu)*SYX                      
        DU( 8)= 1/(2*mu)*(SYY-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                    
        DU( 9)= 1/(2*mu)*SZY                     
        DU(10)= 1/(2*mu)*SZX                     
        DU(11)= 1/(2*mu)*SZY                        
        DU(12)= 1/(2*mu)*(SZZ-lambda/(3*lambda+2*mu)*(SXX+SYY+SZZ))                         
        DO I=1,12                                                   
        U(I)=U(I)+POT3*DU(I) 
        END DO                                       
      ENDIF
        UX=U(1)
        UY=U(2)
        UZ=U(3)
        EXX=U(4)                                                            
        EYX=U(5)                                                            
        EZX=U(6)                                                            
        EXY=U(7)                                                            
        EYY=U(8)                                                            
        EZY=U(9)                                                            
        EXZ=U(10)                                                            
        EYZ=U(11)                                                            
        EZZ=U(12)
      END
