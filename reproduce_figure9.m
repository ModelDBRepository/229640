% Script to reproduce Figure 9
% calls 'closedloop.m' which contains the differential equations   
% initsA lead to tachypnea, initsB lead to eupnea

clear all

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

initsA = [-50.05986089 0.005140176 0.501330626 0.00094653 2.202113749 76.25930796 75.6];
initsB = [-49.69950791 0.005616305 0.528659973 0.000510575 2.126659684 78.26663183 78.1];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[tA,uA] = ode15s('closedloop',[0 420000],initsA,options);

[t1,u1] = ode15s('closedloop',[0 180000],initsB,options);

inits2 = u1(end,:);
inits2(end) = 40;

[t2,u2] = ode15s('closedloop',[180000 360000],inits2,options);

inits3 = u2(end,:);
inits3(end) = 30;

[t3,u3] = ode15s('closedloop',[360000 420000],inits3,options);

tB = [t1; t2; t3];
uB = [u1; u2; u3];

vA=uA(:,1);
vB=uB(:,1);

PO2bloodA = uA(:,7);
PO2bloodB = uB(:,7);

%% Make plot

close('all')

set(0,'DefaultAxesFontSize',24)

lw=2;

figure(1)
plot(tA/1000,tA*0+76.845,'k--','Linewidth',lw)
hold on
plot(tA/1000,PO2bloodA,'r','Linewidth',lw)
plot(tB/1000,PO2bloodB,'b','Linewidth',lw)
xlabel('$t (s)$','Interpreter','latex')
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','latex')
set(gca,'box','off','XTick',0:60:420)
xlim([0 420])





