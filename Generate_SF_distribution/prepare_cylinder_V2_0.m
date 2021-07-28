clear all
%%Massimo Nespoli May 2021
%%% the script generates the distribution of single forces in compatible
%%% format for edcmp_sf.inp
%%% please notice that in this script x is East direction y North
%%%%%%%%%%INPUT%%%%%%%%%%%%
R= 2525;                    % cylinder radius (m)
n= 1000;                    % number of elements on circonference
m= 51;                      % number of elements on depth
q=151;                      %number of elements on each horizontal surface along r (odd number!!!)
tdepth=1425-500/2;          % top depth (m)
bdepth=1425+500/2;          % bottom depth (m)
Vp=2434.3;                  % Vp (m/s)
Vs=1490.7;                  % Vs (m/s)
e0=0.001;                   % potency
rho=2700;                   % Density (kg/m^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=rho*(Vp^2-4/3*Vs^2);
Pressioneass=3*K*e0;

count=0;
for teta=0:360/n:360-360/n
    if m==1
     zeta=(bdepth+tdepth)/2;
        count=count+1;
        number(count)=count;
        tetaQ(count)=teta;
        zetaQ(count)=zeta;
        xQ(count)=R*cosd(teta);
        yQ(count)=R*sind(teta);
        zQ(count)=zeta;
        azimuth(count)=teta;
        inclination(count)=90;
    
    else
        for zeta=tdepth:(bdepth-tdepth)/m:bdepth
        count=count+1;
        number(count)=count;
        tetaQ(count)=teta;
        zetaQ(count)=zeta;
        xQ(count)=R*cosd(teta);
        yQ(count)=R*sind(teta);
        zQ(count)=zeta;
        azimuth(count)=teta;
        inclination(count)=90;
        end
    end
end

%%%%% define pressure for points %%%%%%%%%%%%%
pressure(1:count)=Pressioneass;    % MPa
%%%% the two lasting elements of pressure are upper and lower faces
%%%% pressure
dl=2*pi*R/n;
dh=(bdepth-tdepth)/m;
dS(1:count)=dl*dh;

nsurf=0;
zeta=tdepth;
%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xqline=linspace(min(xQ),max(xQ),q);
yqline=linspace(min(yQ),max(yQ),q);
for t=1:length(xqline)
    for r=1:length(yqline)
        xspace=xqline(t);
        yspace=yqline(r);
        distance=sqrt(xspace^2+yspace^2);
        if distance <=R
                nsurf=nsurf+1;
                count=count+1;
                number(count)=count;
                xQ(count)=xspace;
                yQ(count)=yspace;
                zQ(count)=zeta;
                pressure(count)=Pressioneass;
                if xspace==0 && yspace==0
                    azimuth(count)=1;
                else
                 azimuth(count)=atan(yspace/xspace);
                end
                 inclination(count)=180;
        end
    end
end

zeta=bdepth;
for t=1:length(xqline)
    for r=1:length(yqline)
        xspace=xqline(t);
        yspace=yqline(r);
        distance=sqrt(xspace^2+yspace^2);
        if distance <=R 
                nsurf=nsurf+1;
                count=count+1;
                number(count)=count;
                xQ(count)=xspace;
                yQ(count)=yspace;
                zQ(count)=zeta;
                pressure(count)=Pressioneass;
                if xspace==0 && yspace==0
                    azimuth(count)=1;
                else
                 azimuth(count)=atan(yspace/xspace);
                end
                 inclination(count)=0;
        end
    end
end

dS(end+1:end+nsurf)=pi*((R)^2)/(nsurf/2);

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

matrix=[number' Forces' yQ' xQ' zQ' azimuth'-360 inclination'];
dlmwrite('EDCMP_sf_input.txt',matrix,'delimiter','\t','precision',6)



