% Script to reproduce Figure 12
% calls 'closedloop.m' which contains the differential equations   
%
% For panels A and C (recovery to eupnea) set breakDur = 49246.6
% For panels B and D (failure to tachypnea) set breakDur = 49246.7

clear all

global breakFlag breakVal
 
%% Parameters

% Duration of chemosensory feedback interruption
breakDur = 49246.6; % panels A and C (recovery to eupnea)
%breakDur = 49246.7; % panels B and D (failure to tachypnea)

% Value of gtonic during chemosensory feedback interruption
breakVal = 0.1;

%% Initial Conditions
v0 = -60; n0 = 0; h0 = 0.6; alpha0 = 0; vollung0 = 2; PO2lung0 = 110; PO2blood0 = 110;
inits0=[v0 n0 h0 alpha0 vollung0 PO2lung0 PO2blood0];

tf = 6e4;

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

breakFlag = 0;
[t0,u0] = ode15s('closedloop_clampgtonic',[0 tf],inits0,options);

inits1 = u0(end,:);
[t1,u1] = ode15s('closedloop_clampgtonic',[0 tf],inits1,options);

inits2 = u1(end,:);
breakFlag = 1;

[t2,u2] = ode15s('closedloop_clampgtonic',[tf tf+breakDur],inits2,options);

inits3 = u2(end,:);
breakFlag = 0;

[t3,u3] = ode15s('closedloop_clampgtonic',[tf+breakDur tf+breakDur+3*tf],inits3,options);

t = [t1; t2; t3];
u = [u1; u2; u3];

v = u(:,1);
n =u(:,2);
h =u(:,3);
alpha =u(:,4);
vollung =u(:,5);
PO2lung =u(:,6);
PO2blood =u(:,7);

BCinds=find(t>=0 & t<=tf);
DCinds=find(t>=tf & t<=(tf+breakDur));
ACinds=find(t>=(tf+breakDur));


%% Make plots

close('all')

set(0,'DefaultAxesFontSize',16)

tsec=t/1000; tfsec=max(t)/1000;

xlo=0;
xhi=tfsec;

% panel A or B

figure(1)
plot(tsec(BCinds),PO2blood(BCinds),'k','Linewidth',2)
hold on
plot(tsec(DCinds),PO2blood(DCinds),'b','Linewidth',2)
if PO2blood(end)<80
    plot(tsec(ACinds),PO2blood(ACinds),'r','Linewidth',2)
else
    plot(tsec(ACinds),PO2blood(ACinds),'g','Linewidth',2)
end
xlabel('t (s)')
ylabel('P_aO_2')
xlim([xlo xhi])
set(gca,'box','off','XTick',0:60:240)
grid on

%%

% panel C or D

figure(2)
plot3(h(BCinds),vollung(BCinds),PO2blood(BCinds),'k','Linewidth',2)
hold on
plot3(h(DCinds),vollung(DCinds),PO2blood(DCinds),'b','Linewidth',2)
if PO2blood(end)<80
    plot3(h(ACinds),vollung(ACinds),PO2blood(ACinds),'r','Linewidth',2)
else
    plot3(h(ACinds),vollung(ACinds),PO2blood(ACinds),'g','Linewidth',2)
end
plot3(h(end),vollung(end),PO2blood(end),'ko','MarkerSize',8,'MarkerFaceColor','k')
xlabel('$h$','Interpreter','Latex')
ylabel('vol$_\mathrm{L}$','Interpreter','Latex')
zlabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
grid on



