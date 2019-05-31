% Script to reproduce Figure 1
% calls 'closedloop.m' 

clear all

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

inits = [-56.8172 9.5344e-04 0.7454 2.0026e-04 2.0525 98.9638 97.7927];

tf = 15000;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[t,u] = ode15s('closedloop',[0 tf],inits,options);

v = u(:,1);
alpha = u(:,4);
vollung = u(:,5);
PO2lung = u(:,6);
PO2blood = u(:,7);
gtonic = 0.3*(1-tanh((PO2blood-85)./30));


%% make plots

set(0,'DefaultAxesFontSize',24)

figure(1)

lw=3;

subplot(2,3,1)
plot(t/1000,v,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out')
ylim([-70 20])
ylabel('$V$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,3,2)
plot(t/1000,alpha,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out','YTick',[0:.005:.01],'YTickLabel',[0:.005:.01])
ylim([-0.001 0.01])
ylabel('$\alpha$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,3,3)
plot(t/1000,vollung,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out','YTick',[2:.5:3],'YTickLabel',[2:.5:3])
ylim([1.9 3])
ylabel('$\mathrm{vol}_\mathrm{L}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,3,4)
plot(t/1000,gtonic,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out')
ylabel('$g_\mathrm{tonic}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,3,5)
plot(t/1000,PO2blood,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out')
ylabel('$P_\mathrm{A}\mathrm{O}_2$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

subplot(2,3,6)
plot(t/1000,PO2lung,'k','Linewidth',lw)
xlim([0 15])
set(gca,'box','off','TickDir','out')
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

set(gcf,'position',get(0,'screensize'))


