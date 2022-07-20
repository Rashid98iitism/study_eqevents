%% Intraplate GMPEs 

% GMPEs
% Input parameters
R1=[30:10:300]; R2=[10:10:300]; R3=[10:10:200]; m=6.0;R6=[100:50:1500];
j=1;
%%
% 1. Raghunath and Iyenger (2007)
 for i=1:length(R1)
     % pga_1 in g
    
     pga_1(i)=(exp(1.7236+0.9453*(m-6)-0.0725*(m-6)^2-log(R1(i))-0.0064*R1(i))); 
 end
 
 %%
 % 2. Torro (1997)modified in 2002
 for i=1:length(R2)
     % pga_2 in g
     Rm=sqrt(R2(i)^2+9.3^2);
     Rmmax=max(log(Rm/100),0);
     pga_2(i)=exp(2.20+0.81*(m-6)-1.27*log(Rm)-(1.16-0.0021)*Rmmax-0.0021*Rm);
 end
 
 %%
 % 3. Hwang (1990)
      h=10.0;
 for i=1:length(R6)
     % pga_3 in g
     Rm=0.06*exp(0.7*m);
     pga_3(i)=exp(-2.904+0.926*m-1.271*log(sqrt(R6(i)^2+h^2)+Rm)-0.00302*sqrt(R6(i)^2+h^2));
 end
 
 %% 
 % 4. Dahle et al. (1990)
 for i=1:length(R2)
     
     if R2(i)<=100
        G=1/R2(i);
     else 
        G=0.01*(100/R2(i))^(5/6);
     end
     pga_4(i)=(exp(-1.471+0.849*m-0.00418*R2(i)+log(G)))/9.81;
 end
 
%%
% 5. NDMA (2010)

for ii=1:length(R6)
a=[log((sqrt(R6(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga_5(j)=(exp(-3.7671+1.2303*m-0.0019*(m)^2-0.0027*(sqrt(R6(ii)^2+(2.5)^2))-1.4857*log(sqrt(R6(ii)^2+(2.5)^2)+0.0385*exp(0.8975*m))+0.1301*log10(sqrt(R6(ii)^2+(2.5)^2))*axx));
j=j+1;
end

  
%% 6. MY GMPE (developed for Chhotanagpur Plateau)

% with coefficients of PGA
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_6(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end

%%
plot(R1,pga_1,'b','linewidth',1.5);hold on;
plot(R2,pga_2,'c','linewidth',1.5);hold on;
plot(R6,pga_3,'m','linewidth',1.5); hold on;
plot(R2,pga_4,'g','linewidth',1.5);hold on;
plot(R6,pga_5,'k','linewidth',1.5);hold on;
plot(R6,pga_6,'r','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Raghukanth and Iyenger (2007)','Toro(2002)','Hwang(1990)','Dahle et al.(1990)','NDMA(2010)','Present Study');

%% GMPEs at 0.1 sec period

 %1. Raghunath and Iyenger (2007)
 for i=1:length(R1)
     % pga_1 in g
     pga_1(i)=(exp(2.4597+0.9620*(m-6)-0.0811*(m-6)^2-log(R1(i))-0.0053*R1(i))); 
 end
 
%3. Hwang (1990)
      h=10.0;
 for i=1:length(R6)
     % pga_3 in g
     Rm=0.06*exp(0.7*m);
     pga_3(i)=exp(-2.312+0.924*m-1.233*log(sqrt(R6(i)^2+h^2)+Rm)-0.00317*sqrt(R6(i)^2+h^2));
 end
 
 %5. NDMA (2010)
for ii=1:length(R6)
a=[log((sqrt(R6(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga_5(j)=(exp(-4.2782+1.4909*m-0.0242*(m)^2-0.0026*(sqrt(R6(ii)^2+(2.5)^2))-1.3652*log(sqrt(R6(ii)^2+(2.5)^2)+0.0288*exp(0.9126*m))+0.1140*log10(sqrt(R6(ii)^2+(2.5)^2))*axx));
j=j+1;
end

% 6. MY GMPE (developed for Chhotanagpur Plateau)
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_6(i)=exp(3.5416+0.9811*(m-6)-0.1041*(m-6)^2-log(R(i))-0.0052*R(i)+log(0.3831));
end

plot(R1,pga_1,'b','linewidth',1.5);hold on;
plot(R6,pga_3,'m','linewidth',1.5); hold on;
plot(R6,pga_5,'k','linewidth',1.5);hold on;
plot(R6,pga_6,'r','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Raghukanth and Iyenger (2007)','Hwang(1990)','NDMA(2010)','Present Study');


%% GMPEs at 0.2 sec period

 %1. Raghunath and Iyenger (2007)
 for i=1:length(R1)
     % pga_1 in g
     pga_1(i)=(exp(1.9192+1.0619.*(m-6)-0.1296*(m-6)^2-log(R1(i))-0.0034*R1(i))); 
 end
 
%3. Hwang (1990)
      h=10.0;
 for i=1:length(R6)
     % pga_3 in g
     Rm=0.06*exp(0.7*m);
     pga_3(i)=exp(-2.968+0.952*m-1.219*log(sqrt(R6(i)^2+h^2)+Rm)-0.00240*sqrt(R6(i)^2+h^2));
 end
 
 %5. NDMA (2010)
for ii=1:length(R6)
a=[log((sqrt(R6(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga_5(j)=(exp(-7.5046+2.241*m-0.0829*(m)^2-0.0024*(sqrt(R6(ii)^2+(2.5)^2))-1.3075*log(sqrt(R6(ii)^2+(2.5)^2)+0.0320*exp(0.8878*m))+0.1026*log10(sqrt(R6(ii)^2+(2.5)^2))*axx));
j=j+1;
end

% 6. MY GMPE (developed for Chhotanagpur Plateau)
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_6(i)=exp(3.071+1.0882*(m-6)-0.1546*(m-6)^2-log(R(i))-0.0042*R(i)+log(0.3410));
end

plot(R1,pga_1,'b','linewidth',1.5);hold on;
plot(R6,pga_3,'m','linewidth',1.5); hold on;
plot(R6,pga_5,'k','linewidth',1.5);hold on;
plot(R6,pga_6,'r','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Raghukanth and Iyenger (2007)','Hwang(1990)','NDMA(2010)','Present Study');


%% GMPEs at 0.5 sec period

 %1. Raghunath and Iyenger (2007)
 for i=1:length(R1)
     % pga_1 in g
     pga_1(i)=(exp(1.1638+1.13615*(m-6)-0.2546*(m-6)^2-log(R1(i))-0.0021*R1(i))); 
 end
 
%3. Hwang (1990)
      h=10.0;
 for i=1:length(R6)
     % pga_3 in g
     Rm=0.06*exp(0.7*m);
     pga_3(i)=exp(-4.344+1.048*m-1.180*log(sqrt(R6(i)^2+h^2)+Rm)-0.00200*sqrt(R6(i)^2+h^2));
 end
 
 %5. NDMA (2010)
for ii=1:length(R6)
a=[log((sqrt(R6(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga_5(j)=(exp(-14.6664+4.1214*m-0.2156*(m)^2-0.0022*(sqrt(R6(ii)^2+(2.5)^2))-1.2735*log(sqrt(R6(ii)^2+(2.5)^2)+0.0601*exp(0.7992*m))+0.0987*log10(sqrt(R6(ii)^2+(2.5)^2))*axx));
j=j+1;
end

% 6. MY GMPE (developed for Chhotanagpur Plateau)
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_6(i)=exp(2.3364+1.11245*(m-6)-0.2806*(m-6)^2-log(R(i))-0.0029*R(i)+log(0.3112));
end

plot(R1,pga_1,'b','linewidth',1.5);hold on;
plot(R6,pga_3,'m','linewidth',1.5); hold on;
plot(R6,pga_5,'k','linewidth',1.5);hold on;
plot(R6,pga_6,'r','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Raghukanth and Iyenger (2007)','Hwang(1990)','NDMA(2010)','Present Study');

%% 

plot(R6,pga_6_0,'r','linewidth',1.5);hold on;
plot(R6,pga_6_01,'b','linewidth',1.5);hold on;
plot(R6,pga_6_02,'g','linewidth',1.5);hold on;
plot(R6,pga_6_03,'k','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('PGA','PSA at 0.1g','PSA at 0.2g','PSA at 0.5g');
set(gcf,'color','w');

%% GMPEs at 1 sec period

 %1. Raghunath and Iyenger (2007)
 for i=1:length(R1)
     % pga_1 in g
     pga_1(i)=(exp(0.3604+1.6791*(m-6)-0.3248*(m-6)^2-log(R1(i))-0.0014*R1(i))); 
 end
 
%3. Hwang (1990)
      h=10.0;
 for i=1:length(R6)
     % pga_3 in g
     Rm=0.06*exp(0.7*m);
     pga_3(i)=exp(-6.362+1.231*m-1.1*log(sqrt(R6(i)^2+h^2)+Rm)-0.00243*sqrt(R6(i)^2+h^2));
 end
 
 %5. NDMA (2010)
for ii=1:length(R6)
a=[log((sqrt(R6(ii)^2+(2.5)^2))/100) 0];
axx=max(a);
pga_5(j)=(exp(-20.3915+5.4297*m-0.2996*(m)^2-0.002*(sqrt(R6(ii)^2+(2.5)^2))-1.2441*log(sqrt(R6(ii)^2+(2.5)^2)+0.0512*exp(0.8229*m))+0.0932*log10(sqrt(R6(ii)^2+(2.5)^2))*axx));
j=j+1;
end

% 6. MY GMPE (developed for Chhotanagpur Plateau)
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_6(i)=exp(0.4415+1.3912*(m-6)-0.5095*(m-6)^2-log(R(i))-0.0050*R(i)+log(0.3112));
end

plot(R1,pga_1,'b','linewidth',1.5);hold on;
plot(R6,pga_3,'m','linewidth',1.5); hold on;
plot(R6,pga_5,'k','linewidth',1.5);hold on;
plot(R6,pga_6,'r','linewidth',1.5);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Raghukanth and Iyenger (2007)','Hwang(1990)','NDMA(2010)','Present Study');

%% New GMPE at different magnitudes 

m=3.0;
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_a(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end
% with coefficients of PGA AT M=4.0
m=4.0;
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_b(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end

% with coefficients of PGA
m=5.0;
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_c(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end

% with coefficients of PGA
m=6.0;
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_d(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end

% with coefficients of PGA
m=7.0;
for i=1:length(R6)
    R(i)=sqrt(R6(i)^2+100);
    pga_e(i)=exp(2.6416+0.9631*(m-6)-0.0941*(m-6)^2-log(R(i))-0.0057*R(i)+log(0.4331));
end

plot(R6,pga_a,'r','linewidth',1.5);hold on;
plot(R6,pga_b,'b','linewidth',1.5);hold on;
plot(R6,pga_c,'k','linewidth',1.5); hold on;
plot(R6,pga_d,'m','linewidth',1.5);hold on;
plot(R6,pga_e,'g','linewidth',1.5);hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
%set(gca, 'XAxisLocation', 'top');  % to get x axis on the top 
%set(gca, 'TickDir', 'out');
%set (gca,'Ydir','reverse');
legend('M_w=3.0','M_w=4.0','M_w=5.0','M_w=6.0','M_w=7.0');

%% Crustal Amplification 

