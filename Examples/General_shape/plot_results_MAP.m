clear all
%%Massimo Nespoli May 2021
%%% the script read output files izmhs.disp, izmhs.strss, izmhs.strn
%%% and plots them along x (East) 
%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%
a0=1;                      % Normalization for x plot (m)
FactorAxis_Stress=1e6;        %factor for stress plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


all_data=importdata('izmhs.disp',' ',3);
yedg=all_data.data(:,1);
xedg=all_data.data(:,2);
Uy=all_data.data(:,3);
Ux=all_data.data(:,4);
Uz=all_data.data(:,5);
all_data_stress=importdata('izmhs.strss',' ',3);
sNN=all_data_stress.data(:,3);
sEE=all_data_stress.data(:,4);
sZZ=all_data_stress.data(:,5);
sNE=all_data_stress.data(:,6);
sEZ=all_data_stress.data(:,7);
sZN=all_data_stress.data(:,8);
all_data_stress=importdata('izmhs.strn',' ',3);
eNN=all_data_stress.data(:,3);
eEE=all_data_stress.data(:,4);
eZZ=-all_data_stress.data(:,5);
eNE=all_data_stress.data(:,6);
eEZ=all_data_stress.data(:,7);
eZN=all_data_stress.data(:,8);
xline=xedg(find(yedg==0 & xedg>=0));
yline=yedg(find(yedg==0 & xedg>=0));
eEEline=eEE(find(yedg==0 & xedg>=0));
sEEline=sEE(find(yedg==0 & xedg>=0));
sNNline=sNN(find(yedg==0 & xedg>=0));
sZZline=sZZ(find(yedg==0 & xedg>=0));
eNNline=eNN(find(yedg==0 & xedg>=0));
eZZline=eZZ(find(yedg==0 & xedg>=0));
Uzline=Uz(find(yedg==0 & xedg>=0));
Uxline=Ux(find(yedg==0 & xedg>=0));
Uyline=Uy(find(yedg==0 & xedg>=0));
eNEline=eNE(find(yedg==0 & xedg>=0));
eEZline=eEZ(find(yedg==0 & xedg>=0));
eZNline=eZN(find(yedg==0 & xedg>=0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

plot(xline/a0,(sEEline-(sEEline+sNNline+sZZline)/3),'LineWidth',2);
grid on
hold on
plot(xline/a0,(sZZline-(sEEline+sNNline+sZZline)/3),'LineWidth',2);
plot(xline/a0,(sNNline-(sEEline+sNNline+sZZline)/3),'LineWidth',2);


ylim([-1.5e6 3e6]);
title('Deviatoric stress (Pa)');
xlabel('E/R');
ylabel('Deviatoric stresses (Pa)');
legend('SEE','SNN','szz');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

for ty=1:length(xline)
    %%%%add the cycle if plotting on a plane crossing the TPE inclusion
    % if sqrt(xline(ty)^2+yline(ty)^2)<a0 
    %     sEElinetot(ty)=(sEEline(ty)-(sEEline(ty)+sNNline(ty)+sZZline(ty))/3-4/3*G*e1);
    %     sNNlinetot(ty)=(sNNline(ty)-(sEEline(ty)+sNNline(ty)+sZZline(ty))/3-4/3*G*e1);
    %     sZZlinetot(ty)=(sZZline(ty)-(sEEline(ty)+sNNline(ty)+sZZline(ty))/3-4/3*G*e1);
    % else
        sEElinetot(ty)=(sEEline(ty));
        sNNlinetot(ty)=(sNNline(ty));
        sZZlinetot(ty)=(sZZline(ty));      
    %end
end

plot(xline/a0,sEElinetot,'LineWidth',2);
hold on
grid on
plot(xline/a0,sNNlinetot,'LineWidth',2);
%-2/3*G*e1
plot(xline/a0,sZZlinetot,'LineWidth',2);

title('Total stress (Pa)');
xlabel('E/R');
ylabel('Total stresses (Pa)');
legend('sEE','sNN','szz');

%%%%%%%%%%%%%%%%%%%%%EDGRN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)

uE=Uxline;
uN=Uyline;
uZ=Uzline;
plot(xline/a0,uE,'LineWidth',2);
grid on
hold on
plot(xline/a0,uN,'LineWidth',2);
plot(xline/a0,uZ,'LineWidth',2);

title('Displacement (m)');
xlabel('E/R');
ylabel('Displacement (m)');

legend('uE','uN','uZ');

%%%%%%MAPS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)

xlin2=linspace(min(xedg),max(xedg),500);
ylin2=linspace(min(yedg),max(yedg),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xedg,yedg,Uz,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar

caxis([-0.5 0.5]);
hold on

xlinq=linspace(min(xedg),max(xedg),20);
ylinq=linspace(min(yedg),max(yedg),20);
[X,Y]=meshgrid(xlinq,ylinq);
disp1 =griddata(xedg,yedg,Ux,X,Y); 
disp2 =griddata(xedg,yedg,Uy,X,Y); 

quiver(xlinq,ylinq,disp1,disp2,'color','green' );
title('Displacement (m)');
xlabel('E (m)');
ylabel('N (m)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)

xlin2=linspace(min(xedg),max(xedg),500);
ylin2=linspace(min(yedg),max(yedg),500);
[X,Y]=meshgrid(xlin2,ylin2);
disp =griddata(xedg,yedg,sEE,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar

caxis([-1 1]*FactorAxis_Stress);

title('SEE (Pa)');
xlabel('E (m)');
ylabel('N (m)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)


disp =griddata(xedg,yedg,sNN,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar

caxis([-1 1]*FactorAxis_Stress);

title('SNN (Pa)');
xlabel('E (m)');
ylabel('N (m)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)


disp =griddata(xedg,yedg,sZZ,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar

caxis([-1 1]*FactorAxis_Stress);

title('Szz (Pa)');
xlabel('E (m)');
ylabel('N (m)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)


disp =griddata(xedg,yedg,sEZ,X,Y); 
L=image(xlin2,ylin2,disp,'Cdatamapping','scaled');
 colorbar

caxis([-1 1]*FactorAxis_Stress);

title('SEz (Pa)');
xlabel('E (m)');
ylabel('N (m)');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 save EDCMPsf Ux Uy Uz sNN sEE sZZ sNE sEZ sZN xedg yedg
 FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
 for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(iFig);
  saveas(FigHandle, [FigName,'.jpg']);
end
  
