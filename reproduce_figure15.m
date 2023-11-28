% Script to reproduce Figure 15
% calls 'closedloop.m' which contains the differential equations   
%
% For panels A and C (recovery to eupnea) set breakDur = 24500 or 24597.4
% For panels B and D (failure to tachypnea) set breakDur = 24600 or 24597.5

clear all

global breakFlag breakVal
 
%% Parameters

% Duration of chemosensory feedback interruption
% to get closer to the boundary, use parameter values 24597.4 and 24597.5
breakDur = 24500; % panels A and C (recovery to eupnea)
%breakDur = 24600; % panels B and D (failure to tachypnea)

% Value of gtonic during chemosensory feedback interruption
breakVal = 0.5;

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

set(0,'DefaultAxesFontSize',24)

tsec=t/1000; tfsec=max(t)/1000;

xlo=0;
xhi=tfsec;

% panel A or B

figure(1)
hold on
plot(tsec(DCinds),PO2blood(DCinds),'b','Linewidth',2)
if PO2blood(end)<80
    plot(tsec(ACinds),PO2blood(ACinds),'r','Linewidth',2)
else
    plot(tsec(ACinds),PO2blood(ACinds),'g','Linewidth',2)
end
plot(tsec(BCinds),PO2blood(BCinds),'k','Linewidth',2)
xlabel('t (s)','Interpreter','latex')
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
xlim([xlo xhi])
ylim([20 130])
set(gca,'box','off','XTick',0:60:240,'YTick',20:20:130)
grid on

%%

tachypnea = [0.3481 2.4267 31.0772];

% panel C or D

figure(2)
hold on
plot3(h(DCinds),vollung(DCinds),PO2blood(DCinds),'b','Linewidth',2)
if PO2blood(end)<80
    plot3(h(ACinds),vollung(ACinds),PO2blood(ACinds),'r','Linewidth',2)
else
    plot3(h(ACinds),vollung(ACinds),PO2blood(ACinds),'g','Linewidth',2)
end
plot3(h(BCinds),vollung(BCinds),PO2blood(BCinds),'k','Linewidth',5)
plot3(tachypnea(1),tachypnea(2),tachypnea(3),'ko','MarkerSize',8,'MarkerFaceColor','k')
xlabel('$h$','Interpreter','Latex')
ylabel('vol$_\mathrm{L}$','Interpreter','Latex')
zlabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
grid on
view(3)
xlim([.3 .9])
ylim([2 6])
zlim([20 130])



