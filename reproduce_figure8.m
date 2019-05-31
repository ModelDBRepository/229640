% This script reproduces Figure 8 by simulating the closed-loop and
% open-loop systems over range of M values and finding the average PO2
% blood value

clear all

% Parameters

global M gtonic_open

% Initial Conditions
v0=-60; n0=0; h0=0.6; alpha0=0; vollung0=2; PO2lung0=110; PO2blood0=110;
inits0=[v0 n0 h0 alpha0 vollung0 PO2lung0 PO2blood0];

tf=6e4;

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

M=8e-6;
gtonic_open=0.3;

% closed loop

[t0_closed,u0_closed]=ode15s('closedloopM',[0 tf],inits0,options);

inits1_closed=u0_closed(end,:);
[t1_closed,u1_closed]=ode15s('closedloopM',[0 tf],inits1_closed,options);

% open loop

[t0_open,u0_open]=ode15s('openloopM',[0 tf],inits0,options);

inits1_open=u0_open(end,:);
[t1_open,u1_open]=ode15s('openloopM',[0 tf],inits1_open,options);


Mvals=2e-6:0.1e-6:18e-6;

for ix=1:length(Mvals)

    M=Mvals(ix)
    
    % closed loop
    
    inits2_closed=u1_closed(end,:);
    [t2_closed,u2_closed]=ode15s('closedloopM',[tf 2*tf],inits2_closed,options);

    inits3_closed=u2_closed(end,:);
    [t3_closed,u3_closed]=ode15s('closedloopM',[2*tf 3*tf],inits3_closed,options);

    inits4_closed=u3_closed(end,:);
    [t4_closed,u4_closed]=ode15s('closedloopM',[3*tf 4*tf],inits4_closed,options);

    inits5_closed=u4_closed(end,:);
    [t5_closed,u5_closed]=ode15s('closedloopM',[4*tf 5*tf],inits5_closed,options);

    inits6_closed=u5_closed(end,:);
    [t6_closed,u6_closed]=ode15s('closedloopM',[5*tf 6*tf],inits6_closed,options);

    t_closed=[t1_closed; t2_closed; t3_closed; t4_closed; t5_closed; t6_closed];
    u_closed=[u1_closed; u2_closed; u3_closed; u4_closed; u5_closed; u6_closed];

    po2blood6_closed=u6_closed(:,7);

    avgIntPO2blood_closed(ix)=trapz(t6_closed,po2blood6_closed)/(t6_closed(end)-t6_closed(1));
    
    % open loop

    inits2_open=u1_open(end,:);
    [t2_open,u2_open]=ode15s('openloopM',[tf 2*tf],inits2_open,options);

    inits3_open=u2_open(end,:);
    [t3_open,u3_open]=ode15s('openloopM',[2*tf 3*tf],inits3_open,options);

    inits4_open=u3_open(end,:);
    [t4_open,u4_open]=ode15s('openloopM',[3*tf 4*tf],inits4_open,options);

    inits5_open=u4_open(end,:);
    [t5_open,u5_open]=ode15s('openloopM',[4*tf 5*tf],inits5_open,options);

    inits6_open=u5_open(end,:);
    [t6_open,u6_open]=ode15s('openloopM',[5*tf 6*tf],inits6_open,options);

    t_open=[t1_open; t2_open; t3_open; t4_open; t5_open; t6_open];
    u_open=[u1_open; u2_open; u3_open; u4_open; u5_open; u6_open];

    po2blood6_open=u6_open(:,7);

    avgIntPO2blood_open(ix)=trapz(t6_open,po2blood6_open)/(t6_open(end)-t6_open(1));

end


%% Make plot

set(0,'DefaultAxesFontSize',24)

close(figure(1))

figure(1)
hold on
plot(Mvals,avgIntPO2blood_closed,'k','Linewidth',3)
plot(Mvals,avgIntPO2blood_open,'b','Linewidth',3)
xlim([.2e-5 18e-6])
ylim([1 140])
xlabel('$M$','Interpreter','latex','FontSize',24)
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','latex','Fontsize',24)
h=legend('closed loop','open loop','Location','Northeast');
set(h,'interpreter','latex','FontSize',20)
legend('boxoff')
set(gca,'XTick',[0.4e-5:.4e-5:1.6e-5])
grid on






