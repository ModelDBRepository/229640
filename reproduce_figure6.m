% Script to reproduce Figure 6
% panels A and B call 'closedloop.m' with different initial conditions

clear all

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];
initsB = [-41.7429 0.0313 0.3442 0.0025 2.4355 23.9533 23.3940];

tf = 15000;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[tA,uA] = ode15s('closedloop',[0 tf],initsA,options);
[tB,uB] = ode15s('closedloop',[0 tf],initsB,options);

vA=uA(:,1);
vollungA=uA(:,5);
PO2bloodA=uA(:,7);
gtonicA=0.3*(1-tanh((PO2bloodA-85)./30));

vB=uB(:,1);
vollungB=uB(:,5);
PO2bloodB=uB(:,7);
gtonicB=0.3*(1-tanh((PO2bloodB-85)./30));


%% make plots

set(0,'DefaultAxesFontSize',24)

figure(1)

% eupnea
subplot(4,2,1)
plot(tA/1000,vA,'k','Linewidth',1)
ylim([-65 15])
xlim([0 15])
ylabel('$V$','Interpreter','Latex')
set(gca,'YTick',[-60:30:0])
set(gca,'box','off')

subplot(4,2,3)
plot(tA/1000,vollungA,'k','Linewidth',2)
ylim([1.9 3.1])
xlim([0 15])
ylabel('vol$_\mathrm{L}$','Interpreter','Latex')
set(gca,'box','off')

subplot(4,2,5)
plot(tA/1000,PO2bloodA,'k','Linewidth',2)
ylim([0 120])
xlim([0 15])
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
set(gca,'box','off')

subplot(4,2,7)
plot(tA/1000,gtonicA,'k','Linewidth',2)
ylim([0 0.65])
xlim([0 15])
ylabel('$g_\mathrm{tonic}$','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
set(gca,'YTick',[0:.3:.6])
set(gca,'box','off')

% tachypnea
subplot(4,2,2)
plot(tB/1000,vB,'k','Linewidth',1)
ylim([-65 15])
xlim([0 15])
set(gca,'YTick',[-60:30:0])
set(gca,'box','off')

subplot(4,2,4)
plot(tB/1000,vollungB,'k','Linewidth',2)
ylim([1.9 3.1])
xlim([0 15])
set(gca,'box','off')

subplot(4,2,6)
plot(tB/1000,PO2bloodB,'k','Linewidth',2)
ylim([0 120])
xlim([0 15])
set(gca,'box','off')

subplot(4,2,8)
plot(tB/1000,gtonicB,'k','Linewidth',2)
ylim([0 0.65])
xlim([0 15])
set(gca,'YTick',[0:.3:.6])
xlabel('$t$ (s)','Interpreter','Latex')
set(gca,'box','off')

set(gcf,'position',get(0,'screensize'))


