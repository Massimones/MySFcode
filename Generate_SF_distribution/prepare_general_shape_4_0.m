clear all
%%Massimo Nespoli May 2021
%%% the script generates distributions of single forces for a
%%% non-symmerical 3D geometry
% matlab functions can be found in the official repository
%%% please notice that in this script x is East direction y North
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elsize=50;               % Max dimension of elements (m)
Vp=2434.3;               % Vp (m/s)
Vs=1490.7;               % Vs (m/s)
e0=0.001;                % potency
rho=2700;                % Density (kg/m^3)
refine=0; % 1 if you want to increase the number of forces in each triangle
BPar=0.3; % 0 to 1 for boundaries smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=rho*(Vp^2-4/3*Vs^2);
Pressioneass=3*K*e0;
model = createpde;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%generating an arbitrary 3d shape
rng('default')
xp = (rand([1 100])-0.5)*2000;
yp= (rand([1 100])-0.5)*1000;
zp=(rand([1 100])-0.5)*500+1000;
x=xp;
y=yp;
z=zp;
%  x=xp(sqrt((xp-1000).^2+(yp).^2)<=500 | sqrt((xp+1000).^2+(yp).^2)<=500 |abs(yp)<100); %bow tie
%  y=yp(sqrt((xp-1000).^2+(yp).^2)<=500 | sqrt((xp+1000).^2+(yp).^2)<=500 | abs(yp)<100);%bow tie
%  z=zp(sqrt((xp-1000).^2+(yp).^2)<=500 | sqrt((xp+1000).^2+(yp).^2)<=500 | abs(yp)<100);%bow tie
%   x=xp(sqrt((xp-500).^2+(yp).^2)<=300 | sqrt((xp+500).^2+(yp).^2)<=500 | yp>-500 & yp<0 & abs(xp)<500); %bean
%   y=yp(sqrt((xp-500).^2+(yp).^2)<=300 | sqrt((xp+500).^2+(yp).^2)<=500 | yp>-500 & yp<0 & abs(xp)<500); %bean
%   z=zp(sqrt((xp-500).^2+(yp).^2)<=300 | sqrt((xp+500).^2+(yp).^2)<=500 | yp>-500 & yp<0 & abs(xp)<500); %bean
figure(1000)
scatter3(x,y,z,'fill');
xlabel('E (m)');
ylabel('N (m)');
zlabel('z (m)');
set(gca,'FontSize',18)
set(gca,'fontname','times')
set(gca,'ZDir','reverse')
%%%%%%%%%%%%%%%%%%%%
figure(2000)
DT = delaunayTriangulation(x',y',z');
[T,Xb] = freeBoundary(DT);
BOUN = boundary(x',y',z',BPar);
hh=trisurf(BOUN,x',y',z','Facecolor','red','FaceAlpha',0.1);
xlabel('E (m)');
ylabel('N (m)');
zlabel('z (m)');
set(gca,'FontSize',18)
set(gca,'fontname','times')
set(gca,'ZDir','reverse')

figure(1)
stlwrite('shape.stl',BOUN,[x' y' z'])


geometry=importGeometry(model,'shape.stl');
n_faces=geometry.NumFaces;
mesh = generateMesh(model,'Hmax',elsize,'GeometricOrder','linear');
pdemesh(mesh);
EL=mesh.Elements';
NO=mesh.Nodes';
 axis equal
 hold on  

 xn=NO(:,1);
 yn=NO(:,2);
 zn=NO(:,3);
set(gca,'ZDir','reverse')
 view([-15,15])
 figure(2)
 pdegplot(model,'FaceLabels','on')
 nodes=[];
 for g=1:n_faces 
 nodes{g} = findNodes(mesh,'region','Face',g);    
 end
     
 disp('I am looking for external faces');    
 count=0;
 %%%%%find external faces
 for i=1:length(EL(:,1))
     
     count=count+1;
     ELS=sort(EL(i,:));
     F1=[ELS(:,1), ELS(:,2), ELS(:,3)];
     FACC(count,:)=F1;
     
     for g=1:n_faces
     if sum(nodes{g}==F1(1))>0 && sum(nodes{g}==F1(2))>0 && sum(nodes{g}==F1(3))>0
     LOGIC(count)=1;
     break
     else
     LOGIC(count)=0;
     end
     end
     
     count=count+1;
     F2=[ELS(:,1), ELS(:,2), ELS(:,4)];
     FACC(count,:)=F2;
     
     for g=1:n_faces
     if sum(nodes{g}==F2(1))>0 && sum(nodes{g}==F2(2))>0 && sum(nodes{g}==F2(3))>0
      LOGIC(count)=1;
      break
     else
     LOGIC(count)=0;
     end
     end
     
     count=count+1;
     F3=[ELS(:,2), ELS(:,3), ELS(:,4)];
     FACC(count,:)=F3;
     
     for g=1:n_faces
     if sum(nodes{g}==F3(1))>0 && sum(nodes{g}==F3(2))>0 && sum(nodes{g}==F3(3))>0
    LOGIC(count)=1;
    break
     else
     LOGIC(count)=0;
     end
     end
     
     count=count+1;
     F4=[ELS(:,1), ELS(:,3), ELS(:,4)];
     FACC(count,:)=F4;
     
     for g=1:n_faces
     if sum(nodes{g}==F4(1))>0 && sum(nodes{g}==F4(2))>0 && sum(nodes{g}==F4(3))>0
      LOGIC(count)=1;
      break
     else
     LOGIC(count)=0;
     end
     end
     
     
     
 end
 disp('I am computing the normals'); 
 FACC_app=FACC(logical(LOGIC),:);
 ons=[1 1 1];
 for k=1:length(FACC_app)
    
  V1(k,:)=[NO(FACC_app(k,1),1), NO(FACC_app(k,1),2), NO(FACC_app(k,1),3)];
  V2(k,:)=[NO(FACC_app(k,2),1), NO(FACC_app(k,2),2), NO(FACC_app(k,2),3)];
  V3(k,:)=[NO(FACC_app(k,3),1), NO(FACC_app(k,3),2), NO(FACC_app(k,3),3)];
  normal(k,:)=(cross(V1(k,:)-V2(k,:),V1(k,:)-V3(k,:)))./norm(cross(V1(k,:)-V2(k,:),V1(k,:)-V3(k,:)));
  Baricentro(k,:)=[(V1(k,1)+V2(k,1)+V3(k,1))/3,(V1(k,2)+V2(k,2)+V3(k,2))/3,(V1(k,3)+V2(k,3)+V3(k,3))/3];
  Area(k)=0.5*sqrt(det([[V1(k,1) V2(k,1) V3(k,1)];[V1(k,2) V2(k,2) V3(k,2)];ons])^2 + det([[V1(k,2) V2(k,2) V3(k,2)];[V1(k,3) V2(k,3) V3(k,3)];ons])^2 + det([[V1(k,3) V2(k,3) V3(k,3)];[V1(k,1) V2(k,1) V3(k,1)];ons])^2);
  if refine ==1
  MV1B(k,:)=[(V1(k,1)+Baricentro(k,1))/2,(V1(k,2)+Baricentro(k,2))/2,(V1(k,3)+Baricentro(k,3))/2];
  MV2B(k,:)=[(V2(k,1)+Baricentro(k,1))/2,(V2(k,2)+Baricentro(k,2))/2,(V2(k,3)+Baricentro(k,3))/2];
  MV3B(k,:)=[(V3(k,1)+Baricentro(k,1))/2,(V3(k,2)+Baricentro(k,2))/2,(V3(k,3)+Baricentro(k,3))/2];
  end
 end
 
 figure(3)
  nodifac=cell2mat(nodes);
 NODIB=NO(nodifac,:);

 scatter3(NODIB(:,1),NODIB(:,2),NODIB(:,3),'fill');
  hold on
  
  %%%%%find the normal of each triangle
Bari_object=  [sum(NODIB(:,1))/length(NODIB(:,1)),sum(NODIB(:,2))/length(NODIB(:,2)),sum(NODIB(:,3))/length(NODIB(:,3))];

 %%%%check normal pointing outside
 Sol1=Baricentro+normal*elsize/10;
 Sol2=Baricentro-normal*elsize/10;
 
 INOUT=inpolyhedron(BOUN,[x' y' z'],Sol1);
 INOUT2=inpolyhedron(BOUN,[x' y' z'],Sol2);
 
  
 for i=1:length(normal(:,1))
     if INOUT2(i)==0
         normal(i,:)=-normal(i,:);
     end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter3(Bari_object(:,1),Bari_object(:,2),Bari_object(:,3),'fill');


   
  if refine ==1
   Position=[Baricentro;MV1B;MV2B;MV3B];
   dS=[Area'/4;Area'/4;Area'/4;Area'/4];
   NOR=[normal;normal;normal;normal];
  else
   Position=Baricentro;
   dS=Area;
   NOR=normal;
  end


pressure(1:length(Position(:,1)))=Pressioneass;
Forces=pressure'.*dS;

xforces=Forces.*NOR(:,1);
yforces=Forces.*NOR(:,2);
zforces=Forces.*NOR(:,3);


for i=1:length(zforces)
if xforces(i)==0 && yforces(i)==0
azimuth(i,1)=1;
elseif    xforces(i)<0
azimuth(i,1)=atand(yforces(i)/xforces(i))+180;
else 
azimuth(i,1)=atand(yforces(i)/xforces(i));
end

inclination(i,1)=acosd(zforces(i)/Forces(i)); 
end


xforces2=Forces.*cosd(azimuth).*sind(inclination);
yforces2=Forces.*sind(azimuth).*sind(inclination);
zforces2=Forces.*cosd(inclination);



  figure(4)
  pdemesh(mesh);
  hold on
  quiver3(Position(:,1),Position(:,2),Position(:,3),NOR(:,1),NOR(:,2),NOR(:,3),3);
  set(gca, 'ZDir','reverse')
  number=1:length(xforces);
  matrix=[number' Forces Position(:,2) Position(:,1) Position(:,3) azimuth-360 inclination];
  dlmwrite('EFCMP_input.txt',matrix,'delimiter','\t','precision',6)
     view([-15,15])
     
      FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
 for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = [num2str(iFig),'_s'];
  saveas(FigHandle, [FigName,'.jpg']);
end
     
     
     
 
 

