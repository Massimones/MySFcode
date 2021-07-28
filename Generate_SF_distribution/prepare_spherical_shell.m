clear all

%%Massimo Nespoli May 2021
%%% the script generates the distribution of single forces in compatible
%%% format for edcmp_sf.inp
%%% please notice that in this script x is East direction y North
%%%%%%%%%%INPUT%%%%%%%%%%%%
R= 500;                    % external shell radius (m)
R2= 200;                   % internal shell  radius (m)
n= 150;                    % number of elements on circonference
centerx=0;                 % center x (m)
centery=0;                 % center y (m)
centerz=5000;              % center z (m)
Vp=2434.3;                  % Vp (m/s)
Vs=1490.7;                  % Vs (m/s)
e0=0.001;                   % potency
rho=2700;                   % Density (kg/m^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=rho*(Vp^2-4/3*Vs^2);
Pressioneass=3*K*e0;

count=0;
for teta=0:360/n:360-360/n
    for fi=0+180/(2*n):360/(2*n):180-180/(2*n)
        
        count=count+1;
        inclination(count)=fi;
        azimuth(count)=teta;
       
        xQ(count)=R*cosd(teta)*sind(fi);
        yQ(count)=R*sind(teta)*sind(fi);
        zQ(count)=R*cosd(fi)+centerz;
        zp=R*cosd(fi+360/(4*n))+centerz;
        zm=R*cosd(fi-360/(4*n))+centerz;
        
        dS(count)=abs(2*pi*R*(zp-zm)/n);
        
    end
end
A1=(4*pi*R^2);
A1C=sum(dS);

for teta=0:360/n:360-360/n
    for fi=0+180/(2*n):360/(2*n):180-180/(2*n)
        
        count=count+1;
        inclination(count)=fi;
        azimuth(count)=teta;
       
        xQ(count)=-R2*cosd(teta)*sind(fi);
        yQ(count)=-R2*sind(teta)*sind(fi);
        zQ(count)=-R2*cosd(fi)+centerz;
        
        zp=R2*cosd(fi+360/(4*n))+centerz;
        zm=R2*cosd(fi-360/(4*n))+centerz;
        
        dS(count)=abs(2*pi*R2*(zp-zm)/n);
         
        
    end
end

%%%%% define pressure for points %%%%%%%%%%%%%
pressure(1:count)=Pressioneass;    % MPa
%%%% the two lasting elements of pressure are upper and lower faces
%%%% pressure
A2=(4*pi*R^2+4*pi*R2^2);
A2V=sum(dS);


Forces=pressure.*dS;
xforces=Forces.*cosd(azimuth).*sind(inclination);
yforces=Forces.*sind(azimuth).*sind(inclination);
zforces=Forces.*cosd(inclination);

scatter3(xQ,yQ,zQ,1,'black');
hold on
quiver3(xQ,yQ,zQ,xforces,yforces,zforces);
axis equal

xlabel('E (m)');
ylabel('N (m)');
zlabel('z (m)');

%%%%%%%%%%%%% generate matrix for EDCMP_sf %%%%%%%%%%
% no  Slip  xs ys zs length width azimuth(to east) inclination(depth>0)
number=[1:length(pressure)];
matrix=[number' Forces' yQ' xQ' zQ' azimuth'-360 inclination'];
dlmwrite('EDCMP_sf_input.txt',matrix,'delimiter','\t','precision',6)



