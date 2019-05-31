% Script to reproduce Figure A1
% calls 'closedloop.m'

clear all

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];

tf = 8000;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[t,u] = ode15s('closedloop',[0 tf],initsA,options);

v=u(:,1);
n=u(:,2);
h=u(:,3);
alpha=u(:,4);
vollung=u(:,5);
po2lung=u(:,6);
po2blood=u(:,7);

%% Calculate derivatives

% Parameters

taulb=500;
M=8e-6;

gnap=2.8; gna=28; gk=11.2; gl=2.8; 
Ena=50; Ek=-85; El=-65; Esyn=0;

r=0.001; Tmax=1; VT=2; Kp=5;
E1=0.0025; E2=0.4;

PO2ext=(760-47)*.21; Vol0=2; R=62.364; Temp=310; 
Hb=150; volblood=5; eta=Hb*1.36; gamma=volblood/22400; betaO2=0.03;

c=2.5; K=26;
SaO2=(po2blood.^c)./(po2blood.^c+K^c);
CaO2=eta*SaO2+betaO2*po2blood;
partial=(c*po2blood.^(c-1)).*(1./(po2blood.^c+K^c)-(po2blood.^c)./((po2blood.^c+K^c).^2));

gtonic=0.3*(1-tanh((po2blood-85)./30));

Jlb=(1/taulb)*(po2lung-po2blood).*(vollung/(R*Temp));
Jbt=M*CaO2*gamma;

C=21;  %pF

theta_mp=-40;  %mV
sigma_mp=-6;  %mV
taumax_h=10000;  %ms
theta_h=-48;  %mV
sigma_h=6;  %mV

theta_m=-34;  %mV
sigma_m=-5;  %mV

taumax_n=10;  %ms
theta_n=-29; %mV
sigma_n=-4;  %mV


% Steady state values and time constants
mp_inf=1./(1+exp((v-theta_mp)./sigma_mp));

m_inf=1./(1+exp((v-theta_m)./sigma_m));

h_inf=1./(1+exp((v-theta_h)./sigma_h));
tau_h=taumax_h./cosh((v-theta_h)./(2*sigma_h));

n_inf=1./(1+exp((v-theta_n)./sigma_n));
tau_n=taumax_n./cosh((v-theta_n)./(2*sigma_n));


% Currents

Inap=gnap*mp_inf.*h.*(v-Ena);
Ina=gna*(m_inf.^3).*(1-n).*(v-Ena);
Ik=gk*(n.^4).*(v-Ek);
Il=gl*(v-El);
Itonic=gtonic.*(v-Esyn);

NT=Tmax./(1+exp(-(v-VT)./Kp));

% Lung Volume

dvolrhs=max(0,-E1*(vollung-Vol0)+E2*alpha);

% Differential equations

dv=(-Inap-Ina-Ik-Il-Itonic)/C;
dh=(h_inf-h)./tau_h;
dn=(n_inf-n)./tau_n;
dalpha=r*NT.*(1-alpha)-r*alpha;
dvollung=-E1*(vollung-Vol0)+E2*alpha;
dpo2lung=(1./vollung).*(PO2ext-po2lung).*dvolrhs-Jlb.*(R*Temp./vollung);
dpo2blood=(Jlb-Jbt)./(gamma*(betaO2+eta*partial));

% get max speed

vx_po2blood=max(abs(dpo2blood))/(max(po2blood)-min(po2blood));
vx_vollung=max(abs(dvollung))/(max(vollung)-min(vollung));
vx_po2lung=max(abs(dpo2lung))/(max(po2lung)-min(po2lung));
vx_h=max(abs(dh))/(max(h)-min(h));
vx_alpha=max(abs(dalpha))/(max(alpha)-min(alpha));
vx_v=max(abs(dv))/(max(v)-min(v));
vx_n=max(abs(dn))/(max(n)-min(n));

vx_po2blood_ind=find(abs(dpo2blood)/(max(po2blood)-min(po2blood))==vx_po2blood);
vx_vollung_ind=find(abs(dvollung)/(max(vollung)-min(vollung))==vx_vollung);
vx_po2lung_ind=find(abs(dpo2lung)/(max(po2lung)-min(po2lung))==vx_po2lung);
vx_h_ind=find(abs(dh)/(max(h)-min(h))==vx_h);
vx_alpha_ind=find(abs(dalpha)/(max(alpha)-min(alpha))==vx_alpha);
vx_v_ind=find(abs(dv)/(max(v)-min(v))==vx_v);
vx_n_ind=find(abs(dn)/(max(n)-min(n))==vx_n);


%% make plots

close all

set(0,'DefaultAxesFontSize',24)

lw=3;
ms=10;

figure(1)
subplot(2,3,2)
plot(v,dv/(max(v)-min(v)),'k','Linewidth',lw)
hold on
plot(v(vx_v_ind),vx_v,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$V$','Interpreter','latex')
ylabel('normalized $V''(t)$','Interpreter','latex')
xlim([-65 20])
ylim([-.5 1.2])

subplot(2,3,1)
plot(n,dn/(max(n)-min(n)),'k','Linewidth',lw)
hold on
plot(n(vx_n_ind),vx_n,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$n$','Interpreter','latex')
ylabel('normalized $n''(t)$','Interpreter','latex')
xlim([-.05 1])
ylim([-0.5 2])

subplot(2,3,4)
plot(h,dh/(max(h)-min(h)),'k','Linewidth',lw)
hold on
plot(h(vx_h_ind),-vx_h,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$h$','Interpreter','latex')
ylabel('normalized $h''(t)$','Interpreter','latex')
xlim([0.65 0.78])
ylim([-0.045 0.005])

subplot(2,3,3)
plot(alpha,dalpha/(max(alpha)-min(alpha)),'k','Linewidth',lw)
hold on
plot(alpha(vx_alpha_ind),vx_alpha,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$\alpha$','Interpreter','latex')
ylabel('normalized $\alpha''(t)$','Interpreter','latex')
xlim([-0.001 0.01])
ylim([-.005 .081])

subplot(2,3,5)
plot(vollung,dvollung/(max(vollung)-min(vollung)),'k','Linewidth',lw)
hold on
plot(vollung(vx_vollung_ind),vx_vollung,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$\mathrm{vol}_\mathrm{L}$','Interpreter','latex')
ylabel('normalized $\mathrm{vol}_\mathrm{L}''(t)$','Interpreter','latex')
xlim([1.9 3.1])
ylim([-.0008 .003])

subplot(2,3,6)
plot(po2lung,dpo2lung/(max(po2lung)-min(po2lung)),'b','Linewidth',lw)
hold on
plot(po2blood,dpo2blood/(max(po2blood)-min(po2blood)),'r','Linewidth',lw)
plot(po2lung(vx_po2lung_ind),vx_po2lung,'go','MarkerFaceColor','g','MarkerSize',ms)
plot(po2blood(vx_po2blood_ind),vx_po2blood,'go','MarkerFaceColor','g','MarkerSize',ms)
xlabel('$P\mathrm{O}_2$','Interpreter','latex')
ylabel('normalized $P\mathrm{O}_2''(t)$','Interpreter','latex')
h=legend('$P_\mathrm{A}\mathrm{O}_2$','$P_\mathrm{a}\mathrm{O}_2$','Location','South','Orientation','Horizontal');
legend('boxoff')
set(h,'Interpreter','latex')
ylim([-.0008 .003])
xlim([92 108])

set(gcf,'position',get(0,'screensize'))

