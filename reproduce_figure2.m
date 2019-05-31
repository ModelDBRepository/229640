% Script to reproduce Figure 2
% panel A calls 'closedloop.m' (dynamic gtonic, dynamic h)
% panel B calls 'openloop.m' (static gtonic, dynamic h)
% panel C calls 'closedloop_hstatic.m' (dynamic gtonic, static h)

clear all

global gtonic_open

gtonic_open = 0.3;

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

initsA = [-54.9027 0.0015 0.7519 0.0001 2.0152 95.9033 94.7513];
initsB = [-51.5212 0.0036 0.6120 0.0001 2.0286 85.3022 84.3020];
initsC = [-53.0098 0.0025 0.6000 0.0001 2.0130 86.1800 85.1608];

tf = 8500;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[tA,uA] = ode15s('closedloop',[0 tf],initsA,options);
[tB,uB] = ode15s('openloop',[0 tf],initsB,options);
[tC,uC] = ode15s('closedloop_hstatic',[0 tf],initsC,options);

vA = uA(:,1);
hA = uA(:,3);
PO2bloodA = uA(:,7);
gtonicA = 0.3*(1-tanh((PO2bloodA-85)./30));

vB=uB(:,1);
hB=uB(:,3);
PO2bloodB = uB(:,7);
gtonicB = gtonic_open*ones(length(tB),1);

vC = uC(:,1);
hC = uC(:,3);
PO2bloodC = uC(:,7);
gtonicC = 0.3*(1-tanh((PO2bloodC-85)./30));


%% make plots

set(0,'DefaultAxesFontSize',24)

figure(1)

% panel A
subplot(3,3,1)
plot(tA/1000,vA,'k','Linewidth',3)
xlim([0 7.5])
ylim([-70 20])
set(gca,'box','off','TickDir','out','XTickLabel',[],'XTick',[0:1:8],'YTick',[-60:40:20])
ylabel('$V$','Interpreter','Latex')
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,4)
plot(tA/1000,hA,'k','Linewidth',3)
xlim([0 7.5])
ylim([0.55 0.78])
set(gca,'box','off','TickDir','out','XTickLabel',[],'XTick',[0:1:8],'YTick',[0.55:.1:.75])
ylabel('$h$','Interpreter','Latex')
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,7)
plot(tA/1000,gtonicA,'k','Linewidth',3)
xlim([0 7.5])
ylim([0.1 0.35])
set(gca,'box','off','TickDir','out','XTick',[0:1:8],'YTick',[0.1:.1:.3])
ylabel('$g_\mathrm{tonic}$','Interpreter','Latex')
xlabel('$t$ (s)','Interpreter','Latex')
set(gca,'ticklength',2*get(gca,'ticklength'))


% panel B
subplot(3,3,2)
plot(tB/1000,vB,'Color',[0 0 1],'Linewidth',3)
xlim([0 6])
ylim([-70 20])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTickLabel',[],'XTick',[0:1:8],'YTick',[-60:40:20])
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,5)
plot(tB/1000,hB,'Color',[0 0 1],'Linewidth',3)
xlim([0 6])
ylim([0.55 0.78])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTickLabel',[],'XTick',[0:1:8],'YTick',[0.55:.1:.75])
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,8)
plot(tB/1000,gtonicB,'Color',[0 0 1],'Linewidth',3)
xlim([0 6])
ylim([0.1 0.35])
set(gca,'box','off','TickDir','out','YTickLabel',[],'XTick',[0:1:8])
xlabel('$t$ (s)','Interpreter','Latex')
set(gca,'ticklength',2*get(gca,'ticklength'))

% panel C
subplot(3,3,3)
plot(tC/1000,vC,'Color',[1 0 0],'Linewidth',3)
xlim([0 8.5])
ylim([-70 20])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTickLabel',[],'XTick',[0:1:8],'YTick',[-60:40:20])
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,6)
plot(tC/1000,hC,'Color',[1 0 0],'Linewidth',3)
xlim([0 8.5])
ylim([0.55 0.78])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTickLabel',[],'XTick',[0:1:8],'YTick',[0.55:.1:.75])
set(gca,'ticklength',2*get(gca,'ticklength'))

subplot(3,3,9)
plot(tC/1000,gtonicC,'Color',[1 0 0],'Linewidth',3)
xlim([0 8.5])
ylim([0.1 0.35])
set(gca,'box','off','TickDir','out','YTickLabel',[],'XTick',[0:1:8])
set(gca,'ticklength',2*get(gca,'ticklength'))
xlabel('$t$ (s)','Interpreter','Latex')

set(gcf,'position',get(0,'screensize'))


