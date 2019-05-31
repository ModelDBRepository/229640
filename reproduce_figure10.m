% This script recreates Figure 10 by simulating instantaneous arterial oxygen tension (PaO2) clamp in the closed-loop model, and reading in bifurcation
% curves pre-computed in XPPAUT. To generate the bifurcation curves use the XPP code "generate_figure11_bifurcation_curves.ode"

clear all

%% simulate closed-loop with instantaneous arterial oxygen tension (PaO2) clamp

% Initial conditions [V n h alpha vollung PO2lung PO2blood]
inits1=[-50.001086299999997 0.005214749336000 0.515241494700000 0.000785233164600 2.178562104000000 77.395450929999996 76.9];

t0=0;
tf1=60*1000; tf2=tf1*2; tf3=tf1*3; tf4=tf1*4; tf5=tf1*5; tf6=tf1*6; tf7=tf1*7; tf8=tf1*8; tf9=tf1*9; tf10=tf1*10;

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

[t1,u1]=ode15s('closedloop',[t0 tf1],inits1,options);

inits2=u1(end,:);
[t2,u2]=ode15s('closedloop',[tf1 tf2],inits2,options);

inits3=u2(end,:);
[t3,u3]=ode15s('closedloop',[tf2 tf3],inits3,options);

inits4=u3(end,:);
inits4(end)=40;
[t4,u4]=ode15s('closedloop',[tf3 tf4],inits4,options);

inits5=u4(end,:);
[t5,u5]=ode15s('closedloop',[tf4 tf5],inits5,options);

inits6=u5(end,:);
[t6,u6]=ode15s('closedloop',[tf5 tf6],inits6,options);

inits7=u6(end,:);
inits7(end)=30;
[t7,u7]=ode15s('closedloop',[tf6 tf7],inits7,options);

inits8=u7(end,:);
[t8,u8]=ode15s('closedloop',[tf7 tf8],inits8,options);

inits9=u8(end,:);
[t9,u9]=ode15s('closedloop',[tf8 tf9],inits9,options);

inits10=u9(end,:);
[t10,u10]=ode15s('closedloop',[tf9 tf10],inits10,options);

t=[t1; t2; t3; t4; t5; t6; t7; t8; t9; t10];
u=[u1; u2; u3; u4; u5; u6; u7; u8; u9; u10];

v=u(:,1);
h=u(:,3);
vol=u(:,5);
po2blood=u(:,7);
gtonic=0.3*(1-tanh((po2blood-85)./30));

regburstmin=177270-4900;
regburstmax=178570;

regburstInds=find(t>=regburstmin & t<=regburstmax); 

bigburstmin=t(regburstInds(9679));
bigburstmax=181200;

bigburstInds=find(t>=bigburstmin & t<=bigburstmax);

%% read in pre-computed bifurcation curves (generated in XPPAUT, see generate_figure4_bifurcation_curves.ode)

% panel C top

dataB=dlmread('data_figure4/figure4_panelB_fastsystem_h_bifdiag_gtonic0pt1247.dat'); % cyan
dataC=dlmread('data_figure4/figure4_panelC_fastsystem_h_bifdiag_gtonic0pt2194.dat'); % green
dataD=dlmread('data_figure4/figure4_panelD_fastsystem_h_bifdiag_gtonic0pt1806.dat'); % magenta

% panel C bottom

dataG=dlmread('data_figure10/figure10_panelC_fastsystem_h_bifdiag_gtonic0pt5709.dat'); % green
dataM=dlmread('data_figure10/figure10_panelC_fastsystem_h_bifdiag_gtonic0pt3461.dat'); % magenta

hB=dataB(:,1); hC=dataC(:,1); hD=dataD(:,1); hG=dataG(:,1); hM=dataM(:,1);
vB=dataB(:,2); vC=dataC(:,2); vD=dataD(:,2); vG=dataG(:,2); vM=dataM(:,2);


%% make plots

set(0,'DefaultAxesFontSize',24)

lw=2;

%% panel A

figure(1)

xlo=169000;
xhi=188000;

subplot(4,1,1)
hold on
plot(t/1000,po2blood,'k','Linewidth',lw)
xlim([xlo/1000 xhi/1000])
ylim([30 120])
set(gca,'XTickLabel',[],'YTick',[40 80 120])
ylabel('$P_aO_2$','Interpreter','latex')

subplot(4,1,2)
hold on
plot(t/1000,gtonic,'k','Linewidth',lw)
xlim([xlo/1000 xhi/1000])
ylim([0 .61])
set(gca,'XTickLabel',[],'YTick',[0 .3 .6])
ylabel('$g_\mathrm{tonic}$','Interpreter','latex')

subplot(4,1,3)
hold on
plot(t/1000,v,'k')
xlim([xlo/1000 xhi/1000])
ylim([-70 20])
set(gca,'XTickLabel',[],'YTick',[-60 -20 20])
ylabel('$V$','Interpreter','latex')

subplot(4,1,4)
hold on
plot(t/1000,vol,'k','Linewidth',lw)
xlim([xlo/1000 xhi/1000])
set(gca,'YTick',[2 3.5 5])
ylim([1.9 5.5])
ylabel('$\mathrm{vol}_\mathrm{L}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','Latex')

%% panel B

figure(2)
subplot(2,1,1)
hold on
plot(t/1000,v,'k','Linewidth',lw)
xlim([177.27 178.57])
ylim([-65 15])
ylabel('$V$','Interpreter','latex')
set(gca,'box','off','YTick',[-60:30:0],'XTick',[177.37:.4:178.57],'XTickLabel',{'177.4' '177.8' '178.2' '178.6'})
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,1,2)
hold on
plot(t/1000,v,'k','Linewidth',lw)
xlim([179.9 181.2])
ylim([-65 15])
ylabel('$V$','Interpreter','latex')
set(gca,'box','off','YTick',[-60:30:0],'XTick',[180:.4:181.2],'XTickLabel',{'180.0' '180.4' '180.8' '181.2'})
xlabel('$t$ (s)','Interpreter','latex')

%% panel C

figure(3)
subplot(2,1,1)
plot(hB(2:110),vB(2:110),'c','Linewidth',2)
hold on
plot(hB(110:310),vB(110:310),'k','Linewidth',1)
plot(hC(2:118),vC(2:118),'g','Linewidth',2)
plot(hC(118:318),vC(118:318),'k','Linewidth',1)
plot(hD(2:114),vD(2:114),'m','Linewidth',2)
plot(hD(114:314),vD(114:314),'k','Linewidth',1)
plot(h(regburstInds),v(regburstInds),'k')
xlim([0.2 1])
ylim([-70 10])
xlabel('$h$','Interpreter','latex')
ylabel('$V$','Interpreter','latex')
set(gca,'box','off')

subplot(2,1,2)
plot(hB(2:110),vB(2:110),'c','Linewidth',2)
hold on
plot(hB(110:310),vB(110:310),'k','Linewidth',1)
plot(hG(2:143),vG(2:143),'g','Linewidth',2)
plot(hG(143:243),vG(143:243),'k','Linewidth',1)
plot(hM(2:126),vM(2:126),'m','Linewidth',2)
plot(hM(126:327),vM(126:327),'k','Linewidth',1)
plot(h(bigburstInds),v(bigburstInds),'k')
xlim([0.2 1])
ylim([-70 10])
xlabel('$h$','Interpreter','latex')
ylabel('$V$','Interpreter','latex')
set(gca,'box','off')
