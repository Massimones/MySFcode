clear all
%%Massimo Nespoli May 2021
%%% the script generates the distribution of single forces in compatible
%%% format for edcmp_sf.inp
%%% please notice that in this script x is East direction y North
%%%%%%%%%%INPUT%%%%%%%%%%%%
LE= 2525*2;                 % length E (m)
LN= 2525*2;                 % length N (m)
Num_el_E=200;               % Num of sf along E
Num_el_N=200;               % Num of sf along N
m= 51;                      % number of elements on depth
width=500 ;                 % width (m)
center=1425;                % center depth (m)
Vp=2434.3;                  % Vp (m/s)
Vs=1490.7;                  % Vs (m/s)
e0=0.001;                   % potency
rho=2700;                   % Density (kg/m^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tdepth=center-width/2;                 % top depth (m)
bdepth=center+width/2;               % bottom depth (m)
K=rho*(Vp^2-4/3*Vs^2);
Pressioneass=3*K*e0;
count=0;

%%%upper face
for x=-LE/2:LE/Num_el_E:LE/2
    for y=-LN/2:LN/Num_el_N:LN/2
count=count+1;
number(count)=count;
xQ(count)=x;
yQ(count)=y;
zQ(count)=tdepth;
azimuth(count)=1;
inclination(count)=180;
    end
end


dS(1:count)=LE*LN/count;

%%%bottom face
for x=-LE/2:LE/Num_el_E:LE/2
    for y=-LN/2:LN/Num_el_N:LN/2
count=count+1;
number(count)=count;
xQ(count)=x;
yQ(count)=y;
zQ(count)=bdepth;
azimuth(count)=1;
inclination(count)=0;
    end
end
dS(end+1:count)=dS;

%lateral x left
count2=0;
for z=tdepth:(bdepth-tdepth)/m:bdepth
    for y=-LN/2:LN/Num_el_N:LN/2
        count2=count2+1;
count=count+1;
number(count)=count;
xQ(count)=-LE/2;
yQ(count)=y;
zQ(count)=z;
azimuth(count)=180;
inclination(count)=90;
    end
end

dS(end+1:end+count2)=LN*width/count2;

%lateral x right
for z=tdepth:(bdepth-tdepth)/m:bdepth
    for y=-LN/2:LN/Num_el_N:LN/2
count=count+1;
number(count)=count;
xQ(count)=LE/2;
yQ(count)=y;
zQ(count)=z;
azimuth(count)=0;
inclination(count)=90;
    end
end

dS(end+1:end+count2)=LN*width/count2;

%lateral near
count3=0;
for x=-LE/2:LE/Num_el_E:LE/2
    for z=tdepth:(bdepth-tdepth)/m:bdepth
        count3=count3+1;
count=count+1;
number(count)=count;
xQ(count)=x;
yQ(count)=-LN/2;
zQ(count)=z;
azimuth(count)=-90;
inclination(count)=90;
    end
end
dS(end+1:end+count3)=LE*width/count3;


%lateral far
for x=-LE/2:LE/Num_el_E:LE/2
    for z=tdepth:(bdepth-tdepth)/m:bdepth
count=count+1;
number(count)=count;
xQ(count)=x;
yQ(count)=LN/2;
zQ(count)=z;
azimuth(count)=90;
inclination(count)=90;
    end
end
dS(end+1:end+count3)=LE*width/count3;

%%%%% define pressure for points %%%%%%%%%%%%%
pressure(1:count)=Pressioneass;    % MPa
%%%% the two lasting elements of pressure are upper and lower faces
%%%% pressure

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



