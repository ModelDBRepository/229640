% This script recreates Figure 4 by simulating open- and closed-loop trajectories and reading in bifurcation
% curves pre-computed in XPPAUT. To generate the bifurcation curves use the XPP code "generate_figure4_bifurcation_curves.ode"

clear all

%% simulate open-loop trajectory for panel A

global gtonic_open

gtonic_open = 0.3;

% Initial conditions [V n h alpha vollung PO2lung PO2blood]
inits1=[-55 0.001 0.74 0.0001 2 100 100];

t0=0;
tf1=60*1000; tf2=tf1*2; tf3=tf1*3; tf4=tf1*4; tf5=tf1*5; tf6=tf1*6;

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

[t1,u1]=ode15s('openloop',[t0 tf1],inits1,options);

inits2=u1(end,:);

[t2,u2]=ode15s('openloop',[tf1 tf2],inits2,options);

inits3=u2(end,:);

[t3,u3]=ode15s('openloop',[tf2 tf3],inits3,options);

inits4=u3(end,:);

[t4,u4]=ode15s('openloop',[tf3 tf4],inits4,options);

inits5=u4(end,:);

[t5,u5]=ode15s('openloop',[tf4 tf5],inits5,options);

inits6=u5(end,:);

[t6,u6]=ode15s('openloop',[tf5 tf6],inits6,options);

vopen=u6(:,1);
hopen=u6(:,3);


%% simulate closed-loop trajectory for panels B-D

% Initial conditions [V n h alpha vollung PO2lung PO2blood]
inits1=[-50.001086299999997 0.005214749336000 0.515241494700000 0.000785233164600 2.178562104000000 77.395450929999996 76.9];

tf7=tf1*7; tf8=tf1*8; tf9=tf1*9; tf10=tf1*10;

[t1,u1]=ode15s('closedloop',[t0 tf1],inits1,options);

inits2=u1(end,:);
[t2,u2]=ode15s('closedloop',[tf1 tf2],inits2,options);

inits3=u2(end,:);
[t3,u3]=ode15s('closedloop',[tf2 tf3],inits3,options);

inits4=u3(end,:);
[t4,u4]=ode15s('closedloop',[tf3 tf4],inits4,options);

inits5=u4(end,:);
[t5,u5]=ode15s('closedloop',[tf4 tf5],inits5,options);

inits6=u5(end,:);
[t6,u6]=ode15s('closedloop',[tf5 tf6],inits6,options);

inits7=u6(end,:);
[t7,u7]=ode15s('closedloop',[tf6 tf7],inits7,options);

inits8=u7(end,:);
[t8,u8]=ode15s('closedloop',[tf7 tf8],inits8,options);

inits9=u8(end,:);
[t9,u9]=ode15s('closedloop',[tf8 tf9],inits9,options);

inits10=u9(end,:);
[t10,u10]=ode15s('closedloop',[tf9 tf10],inits10,options);

tclosed=[t1; t2; t3; t4; t5; t6; t7; t8; t9; t10];
uclosed=[u1; u2; u3; u4; u5; u6; u7; u8; u9; u10];

vclosed=uclosed(:,1);
hclosed=uclosed(:,3);

regburstmin=177270-4900;
regburstmax=178570;

regburstInds=find(tclosed>=regburstmin & tclosed<=regburstmax); 

%% read in pre-computed bifurcation curves (generated in XPPAUT, see generate_figure4_bifurcation_curves.ode)

dataA=dlmread('data_figure4/figure4_panelA_fastsystem_h_bifdiag_gtonic0pt3.dat'); % panel A 
dataB=dlmread('data_figure4/figure4_panelB_fastsystem_h_bifdiag_gtonic0pt1247.dat'); % panel B
dataC=dlmread('data_figure4/figure4_panelC_fastsystem_h_bifdiag_gtonic0pt2194.dat'); % panel C
dataD=dlmread('data_figure4/figure4_panelD_fastsystem_h_bifdiag_gtonic0pt1806.dat'); % panel D

hA=dataA(:,1); hB=dataB(:,1); hC=dataC(:,1); hD=dataD(:,1);
vA=dataA(:,2); vB=dataB(:,2); vC=dataC(:,2); vD=dataD(:,2);

% calculate h-nullcline
vrange=-70:.1:0;
theta_h=-48; sigma_h=6; 
hnull=1./(1+exp((vrange-theta_h)./sigma_h));


%% make plot

set(0,'DefaultAxesFontSize',24)

figure(1)

% panel A
subplot(2,4,1)
hold on
plot(hA(2:124),vA(2:124),'k','Linewidth',3)
plot(hA(124:417),vA(124:417),'k','Linewidth',1)
plot(hA(417:end),vA(417:end),'k','Linewidth',3)
plot(hopen,vopen,'b','Linewidth',1)
xlim([0.53 .93])
ylim([-70 10])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'XTick',[.53 .73 .93])
title('Open loop ($g_\mathrm{tonic}$ = 0.3)','FontSize',18,'Interpreter','Latex')

subplot(2,4,5)
hold on
plot(hA(2:124),vA(2:124),'k','Linewidth',3)
plot(hA(124:324),vA(124:324),'k','Linewidth',1)
plot(hopen,vopen,'b','Linewidth',1)
plot(hnull,vrange,'--','Color',[.5 .5 .5],'Linewidth',3)
xlim([0.56 .62])
ylim([-55 -47])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
plot(.6090,-50.6576,'o','Color',[.5 .5 .5],'MarkerSize',10,'Linewidth',3)
set(gca,'XTick',[.56 .59 .62])

% panel B
subplot(2,4,2)
plot(hB(2:110),vB(2:110),'c','Linewidth',3)
hold on
plot(hB(110:444),vB(110:444),'k','Linewidth',1)
plot(hB(444:end),vB(444:end),'c','Linewidth',3)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(regburstInds(9679)),vclosed(regburstInds(9679)),'co','MarkerFaceColor','c','MarkerSize',10)
xlim([0.6 1])
ylim([-70 10])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off')
title('Closed loop ($g_\mathrm{tonic}$ = 0.12)','FontSize',18,'Interpreter','Latex')

subplot(2,4,6)
plot(hB(2:110),vB(2:110),'c','Linewidth',3)
hold on
plot(hB(110:310),vB(110:310),'k','Linewidth',1)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(regburstInds(9679)),vclosed(regburstInds(9679)),'co','MarkerFaceColor','c','MarkerSize',10)
xlim([0.65 0.9])
ylim([-62 -45])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off','XTick',[.65 .75 .85])
plot(hnull,vrange,'--','Color',[.5 .5 .5],'Linewidth',3)
plot(.8457,-58.2080,'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

% panel C
subplot(2,4,3)
hold on
plot(hC(2:118),vC(2:118),'g','Linewidth',3)
plot(hC(118:430),vC(118:430),'k','Linewidth',1)
plot(hC(430:end),vC(430:end),'g','Linewidth',3)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(292457),vclosed(292457),'go','MarkerFaceColor','g','MarkerSize',10)
xlim([0.6 1])
ylim([-70 10])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off','XTick',[.6 .8 1])
title('Closed loop ($g_\mathrm{tonic}$ = 0.22)','FontSize',18,'Interpreter','Latex')

subplot(2,4,7)
hold on
plot(hC(2:118),vC(2:118),'g','Linewidth',3)
plot(hC(118:318),vC(118:318),'k','Linewidth',1)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(292457),vclosed(292457),'go','MarkerFaceColor','g','MarkerSize',10)
xlim([0.65 0.9])
ylim([-62 -45])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off','XTick',[.65 .75 .85])
plot(hnull,vrange,'--','Color',[.5 .5 .5],'Linewidth',3)
plot(.7343,-54.0988,'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

% panel D
subplot(2,4,4)
hold on
plot(hD(2:114),vD(2:114),'m','Linewidth',3)
plot(hD(114:435),vD(114:435),'k','Linewidth',1)
plot(hD(435:end),vD(435:end),'m','Linewidth',3)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(301655),vclosed(301655),'mo','MarkerFaceColor','m','MarkerSize',10)
xlim([0.6 1])
ylim([-70 10])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off')
title('Closed loop ($g_\mathrm{tonic}$ = 0.18)','FontSize',18,'Interpreter','Latex')

subplot(2,4,8)
hold on
plot(hD(2:114),vD(2:114),'m','Linewidth',3)
plot(hD(114:314),vD(114:314),'k','Linewidth',1)
plot(hclosed(regburstInds),vclosed(regburstInds),'k')
plot(hclosed(301655),vclosed(301655),'mo','MarkerFaceColor','m','MarkerSize',10)
xlim([0.65 0.9])
ylim([-62 -45])
xlabel('$h$','Interpreter','Latex')
ylabel('$V$','Interpreter','Latex')
set(gca,'box','off','XTick',[.65 .75 .85])
plot(hnull,vrange,'--','Color',[.5 .5 .5],'Linewidth',3)
plot(.7875,-55.8586,'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)
