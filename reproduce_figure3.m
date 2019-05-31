% This script reproduces Figure 3 by simulating the open-loop system for a 
% range of gtonic values and the closed-loop static h system for a range of
% h values. For the open-loop simulations it finds the minimum and maximum h
% values observed and for the closed-loop static h simulations it finds the minimum
% and maximum gtonic values observed. It also plots the boundaries of
% quiscence, bursting, and beating regimes. These boundaries can be found
% by inspecting the spike time statistics of the simulations.

clear all; close all

tf=12e4;

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

%% open-loop

global gtonic_open

gtonic_vals=0.05:0.01:0.7;

% Initial conditions [V n h alpha vollung PO2lung PO2blood]
inits0 = [-60 0 0 0 2 110 110];

for gx=1:length(gtonic_vals)
    
    gtonic_open = gtonic_vals(gx)

    [t0,u0]=ode15s('openloop',[0 tf],inits0,options);

    inits1=u0(end,:);
    [t1,u1]=ode15s('openloop',[0 tf],inits1,options);
    
    t=t1;
    v=u1(:,1);
    h=u1(:,3);
    
    % find spike time statistics
      
    spikeTimes=[];
    thresh=0;

    for jx=2:length(t)
        if v(jx)>=thresh
            if v(jx-1)<thresh
                spikeTimes=[spikeTimes t(jx)];
            end
        end
    end
    
    numSpikes_open(gx)=length(spikeTimes);
    
    if length(spikeTimes)>1
    
        isi=diff(spikeTimes);
        maxisi_open(gx)=max(isi);
        cv_open(gx)=std(isi)/mean(isi);
    
    else
        
        maxisi_open(gx)=NaN;
        cv_open(gx)=NaN;
        
    end
    
    % find range of h values traversed
    
    hrange_open(gx,:)=[min(h) max(h)];
    
end

%% closed loop static h

hvals=0.2:0.01:0.9;

for hx=1:length(hvals)
    
    h0=hvals(hx)
    
    % Initial conditions [V n h alpha vollung PO2lung PO2blood]
    inits0 = [-60 0 h0 0 2 110 110];
    
    [t0,u0]=ode15s('closedloop_hstatic',[0 tf],inits0,options);

    inits1=u0(end,:);
    [t1,u1]=ode15s('closedloop_hstatic',[0 tf],inits1,options);

    inits2=u1(end,:);
    [t2,u2]=ode15s('closedloop_hstatic',[0 15000],inits2,options);
    
    t=t2;
    v=u2(:,1);
    PO2blood=u2(:,7);
    gtonic=0.3*(1-tanh((PO2blood-85)./30));
    
    % find spike time statistics
      
    spikeTimes=[];
    thresh=0;

    for jx=2:length(t)
        if v(jx)>=thresh
            if v(jx-1)<thresh
                spikeTimes=[spikeTimes t(jx)];
            end
        end
    end
    
    numSpikes_hstatic(hx)=length(spikeTimes);
    
    if length(spikeTimes)>1
    
        isi=diff(spikeTimes);
        maxisi_hstatic(hx)=max(isi);
        cv_hstatic(hx)=std(isi)/mean(isi);
        
    else
        
        maxisi_hstatic(hx)=NaN;
        cv_hstatic(hx)=NaN;
        
    end
    
    % find range of gtonic values traversed
    
    grange_hstatic(hx,:)=[min(gtonic) max(gtonic)];
   
end

%% Initial conditions [V n h alpha vollung PO2lung PO2blood]

initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];
initsB = [-41.7429 0.0313 0.3442 0.0025 2.4355 23.9533 23.3940];

tf = 15000;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

[tA,uA] = ode15s('closedloop',[0 tf],initsA,options);
[tB,uB] = ode15s('closedloop',[0 tf],initsB,options);

hA=uA(:,3);
PO2bloodA=uA(:,7);
gtonicA=0.3*(1-tanh((PO2bloodA-85)./30));

hB=uB(:,3);
PO2bloodB=uB(:,7);
gtonicB=0.3*(1-tanh((PO2bloodB-85)./30));


%% make plot

figure(1)

hold on
% boundaries of Q/Be/Bu/Be for closed-loop static h
plot([0 1],[.35 .35],'r--','Linewidth',1)
plot([0 1],[.46 .46],'r--','Linewidth',1)
plot([0 1],[.78 .78],'r--','Linewidth',1)

% boundaries of Q/Bu/Be for open-loop
plot([.27 .27],[0 1],'b--','Linewidth',1)
plot([.46 .46],[0 1],'b--','Linewidth',1)

% min and max h values during open-loop
plot(gtonic_vals,hrange_open(:,1),'b','Linewidth',3)
plot(gtonic_vals,hrange_open(:,2),'b','Linewidth',3)

% min and max gtonic values during closed-loop static h
plot(grange_hstatic(:,1),hvals,'Color',[1 0 0],'Linewidth',3)
plot(grange_hstatic(:,2),hvals,'Color',[1 0 0],'Linewidth',3)

% closed-loop (dynamic h) eupneic limit cycle
plot(gtonicA,hA,'k','Linewidth',5)

% closed-loop (dynamic h) tachypneic fixed point
plot(gtonicB,hB,'ko','MarkerSize',10,'MarkerFaceColor','k')

gtonicA_max_ind=find(gtonicA==max(gtonicA));
gtonicA_min_ind=find(gtonicA==min(gtonicA));
gtonicA_0pt18_ind=find(abs(gtonicA-0.18)==min(abs(gtonicA-0.18)));

plot(gtonicA(gtonicA_max_ind),hA(gtonicA_max_ind),'go','MarkerSize',10,'MarkerFaceColor','g')
plot(gtonicA(gtonicA_min_ind),hA(gtonicA_min_ind),'co','MarkerSize',10,'MarkerFaceColor','c')
plot(gtonicA(gtonicA_0pt18_ind),hA(gtonicA_0pt18_ind),'mo','MarkerSize',10,'MarkerFaceColor','m')

ylim([0.2 0.9])
xlim([0.05 0.7])

ylabel('$h$','Interpreter','Latex')
xlabel('$g_\mathrm{tonic}$','Interpreter','Latex')
set(gca,'YTick',[0.1:.1:1])




    
    

    