% This script reproduces Figure 5
% first the open-loop, then closed-loop, then gtonic forcing

clear all

scale_factors=[0.8 1 1.5 2 2.5 3 3.5 4 4.5];

%% initialize open and closed-loop

global taumax_h gtonic_open

v0=-55; n0=0.001; h0=0.74; alpha0=0.0001; vollung0=2; PO2lung0=100; PO2blood0=100; 

inits1_open_closed=[v0 n0 h0 alpha0 vollung0 PO2lung0 PO2blood0];

t0=0;
tf1=60*1000;
tf2=60*1000*2;
tf3=60*1000*3;
tf4=60*1000*4;
tf5=60*1000*5;
tf6=60*1000*6;

options=odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@events_figure5);

gtonic_open=0.3;


%% initialize forced gtonic

inits1_forced=[-56.031664239999998 0.001160025770000 0.753151682700000 0.000095924685630 2.024680252000000 97.636469689999998 96.467834690000004];

tf=(6.372270214e3)*20;
tin=0:0.01:tf;

[t1,u1,te1,ye1,ie1]=ode15s(@closedloop,[0 tf],inits1_forced,options);

data=[t1 u1];


%% loop through scale factors

for scalex=1:length(scale_factors)
    
    scale_factor=scale_factors(scalex)
       
    taumax_h=scale_factor*10000;

    % open-loop

    [t1_open,u1_open,te1_open,ye1_open,ie1_open]=ode15s(@openloop_tauh,[t0 tf1],inits1_open_closed,options);

    inits2_open=u1_open(end,:);
    [t2_open,u2_open,te2_open,ye2_open,ie2_open]=ode15s(@openloop_tauh,[t0 tf2],inits2_open,options);

    inits3_open=u2_open(end,:);
    [t3_open,u3_open,te3_open,ye3_open,ie3_open]=ode15s(@openloop_tauh,[t0 tf3],inits3_open,options);

    inits4_open=u3_open(end,:);
    [t4_open,u4_open,te4_open,ye4_open,ie4_open]=ode15s(@openloop_tauh,[t0 tf4],inits4_open,options);

    inits5_open=u4_open(end,:);
    [t5_open,u5_open,te5_open,ye5_open,ie5_open]=ode15s(@openloop_tauh,[t0 tf5],inits5_open,options);

    inits6_open=u5_open(end,:);
    [t6_open,u6_open,te6_open,ye6_open,ie6_open]=ode15s(@openloop_tauh,[t0 tf6],inits6_open,options);

    % closed-loop

    [t1_closed,u1_closed,te1_closed,ye1_closed,ie1_closed]=ode15s(@closedloop_tauh,[t0 tf1],inits1_open_closed,options);

    inits2_closed=u1_closed(end,:);
    [t2_closed,u2_closed,te2_closed,ye2_closed,ie2_closed]=ode15s(@closedloop_tauh,[t0 tf2],inits2_closed,options);

    inits3_closed=u2_closed(end,:);
    [t3_closed,u3_closed,te3_closed,ye3_closed,ie3_closed]=ode15s(@closedloop_tauh,[t0 tf3],inits3_closed,options);

    inits4_closed=u3_closed(end,:);
    [t4_closed,u4_closed,te4_closed,ye4_closed,ie4_closed]=ode15s(@closedloop_tauh,[t0 tf4],inits4_closed,options);

    inits5_closed=u4_closed(end,:);
    [t5_closed,u5_closed,te5_closed,ye5_closed,ie5_closed]=ode15s(@closedloop_tauh,[t0 tf5],inits5_closed,options);

    inits6_closed=u5_closed(end,:);
    [t6_closed,u6_closed,te6_closed,ye6_closed,ie6_closed]=ode15s(@closedloop_tauh,[t0 tf6],inits6_closed,options);

    % forced gtonic
    
    taumax_h=10000;

    tdata=scale_factor*data(:,1);
    vdata=data(:,2);
    ndata=data(:,3);
    hdata=data(:,4);
    po2blood=data(:,end);
    gtonic_orig=0.3*(1-tanh((po2blood-85)./30));

    % parameters
    gnap=2.8; gna=28; gk=11.2; gl=2.8; 
    Ena=50; Ek=-85; El=-65; Esyn=0;

    C=21;  %pF

    theta_mp=-40;  %mV
    sigma_mp=-6;  %mV
    theta_h=-48;  %mV
    sigma_h=6;  %mV

    theta_m=-34;  %mV
    sigma_m=-5;  %mV

    taumax_n=10;  %ms
    theta_n=-29; %mV
    sigma_n=-4;  %mV

    mp_inf=@(v) 1/(1+exp((v-theta_mp)/sigma_mp));
    m_inf=@(v) 1/(1+exp((v-theta_m)/sigma_m));
    h_inf=@(v) 1/(1+exp((v-theta_h)/sigma_h));
    tau_h=@(v) taumax_h/cosh((v-theta_h)/(2*sigma_h));
    n_inf=@(v) 1/(1+exp((v-theta_n)/sigma_n));
    tau_n=@(v) taumax_n/cosh((v-theta_n)/(2*sigma_n));

    tf=max(tdata);
    dt=0.01;

    t=0:dt:tf;

    gtonic=interp1(tdata,gtonic_orig,t);

    numSteps=length(t)-1;

    v=zeros(length(t),1);
    n=zeros(length(t),1);
    h=zeros(length(t),1);

    v(1)=vdata(1);
    n(1)=ndata(1);
    h(1)=hdata(1);

    for ix=1:numSteps

        Inap=gnap*mp_inf(v(ix))*h(ix)*(v(ix)-Ena);
        Ina=gna*(m_inf(v(ix))^3)*(1-n(ix))*(v(ix)-Ena);
        Ik=gk*(n(ix)^4)*(v(ix)-Ek);
        Il=gl*(v(ix)-El);
        Itonic=gtonic(ix)*(v(ix)-Esyn);

        k1v=(-Inap-Ina-Ik-Il-Itonic)/C;
        k1n=(n_inf(v(ix))-n(ix))/tau_n(v(ix));
        k1h=(h_inf(v(ix))-h(ix))/tau_h(v(ix));

        av=v(ix)+k1v*dt;
        an=n(ix)+k1n*dt;
        ah=h(ix)+k1h*dt;

        Inap=gnap*mp_inf(av)*ah*(av-Ena);
        Ina=gna*(m_inf(av)^3)*(1-an)*(av-Ena);
        Ik=gk*(an^4)*(av-Ek);
        Il=gl*(av-El);
        Itonic=gtonic(ix+1)*(av-Esyn);

        k2v=(-Inap-Ina-Ik-Il-Itonic)/C;
        k2n=(n_inf(av)-an)/tau_n(av);
        k2h=(h_inf(av)-ah)/tau_h(av);

        v(ix+1)=v(ix)+(k1v+k2v)*dt/2;
        n(ix+1)=n(ix)+(k1n+k2n)*dt/2;
        h(ix+1)=h(ix)+(k1h+k2h)*dt/2;

    end
    
    t_forced=t;
    v_forced=v;


    %%% find spike times and burst properties

    %% open-loop

    t=t6_open;
    v=u6_open(:,1);

    numT=length(t);
    spikeTimes=[];
    tlap=0;

    for ix=2:numT
        if v(ix)>0
            if v(ix-1)<0
                if (t(ix)-tlap)>10
                    spikeTimes=[spikeTimes t(ix)];
                    tlap=t(ix);
                end
            end
        end
    end

    numSpikes=length(spikeTimes);
    isi=spikeTimes(2:end)-spikeTimes(1:(end-1));

    max_isi=max(isi);

    ibi_inds=find(isi>0.8*max_isi);
    ibi=isi(ibi_inds);
    num_ibi=length(ibi);
    mean_ibi=mean(ibi/1000);

    first_spike_inds=ibi_inds+1;
    first_spikeTimes=spikeTimes(first_spike_inds);

    last_spike_inds=ibi_inds;
    last_spikeTimes=spikeTimes(last_spike_inds);

    period_first=first_spikeTimes(2:end)-first_spikeTimes(1:(end-1));
    period_last=last_spikeTimes(2:end)-last_spikeTimes(1:(end-1));

    mean_period=mean([period_first period_last])/1000;
    mean_burst_dur=mean_period-mean_ibi;

    spikes_per_burst=length(first_spike_inds(end-1):last_spike_inds(end));

    openIBI(scalex)=mean_ibi;
    openPeriod(scalex)=mean_period;
    openBurstDur(scalex)=mean_burst_dur;
    openSPB(scalex)=spikes_per_burst;

    %% closed-loop

    t=t6_closed;
    v=u6_closed(:,1);

    numT=length(t);
    spikeTimes=[];
    tlap=0;

    for ix=2:numT
        if v(ix)>0
            if v(ix-1)<0
                if (t(ix)-tlap)>10
                    spikeTimes=[spikeTimes t(ix)];
                    tlap=t(ix);
                end
            end
        end
    end

    numSpikes=length(spikeTimes);
    isi=spikeTimes(2:end)-spikeTimes(1:(end-1));

    max_isi=max(isi);

    ibi_inds=find(isi>0.8*max_isi);
    ibi=isi(ibi_inds);
    num_ibi=length(ibi);
    mean_ibi=mean(ibi/1000);

    first_spike_inds=ibi_inds+1;
    first_spikeTimes=spikeTimes(first_spike_inds);

    last_spike_inds=ibi_inds;
    last_spikeTimes=spikeTimes(last_spike_inds);

    period_first=first_spikeTimes(2:end)-first_spikeTimes(1:(end-1));
    period_last=last_spikeTimes(2:end)-last_spikeTimes(1:(end-1));

    mean_period=mean([period_first period_last])/1000;
    mean_burst_dur=mean_period-mean_ibi;

    spikes_per_burst=length(first_spike_inds(end-1):last_spike_inds(end));

    closedIBI(scalex)=mean_ibi;
    closedPeriod(scalex)=mean_period;
    closedBurstDur(scalex)=mean_burst_dur;
    closedSPB(scalex)=spikes_per_burst;


    %% gtonic forcing
    
    t=t_forced;
    v=v_forced;

    numT=length(t);
    spikeTimes=[];
    tlap=0;

    for ix=2:numT
        if v(ix)>0
            if v(ix-1)<0
                if (t(ix)-tlap)>10
                    spikeTimes=[spikeTimes t(ix)];
                    tlap=t(ix);
                end
            end
        end
    end

    numSpikes=length(spikeTimes);
    isi=spikeTimes(2:end)-spikeTimes(1:(end-1));

    max_isi=max(isi);

    ibi_inds=find(isi>0.8*max_isi);
    ibi=isi(ibi_inds);
    num_ibi=length(ibi);
    mean_ibi=mean(ibi/1000);

    first_spike_inds=ibi_inds+1;
    first_spikeTimes=spikeTimes(first_spike_inds);

    last_spike_inds=ibi_inds;
    last_spikeTimes=spikeTimes(last_spike_inds);

    period_first=first_spikeTimes(2:end)-first_spikeTimes(1:(end-1));
    period_last=last_spikeTimes(2:end)-last_spikeTimes(1:(end-1));

    mean_period=mean([period_first period_last])/1000;
    mean_burst_dur=mean_period-mean_ibi;

    spikes_per_burst=length(first_spike_inds(end-1):last_spike_inds(end));

    forcedIBI(scalex)=mean_ibi;
    forcedPeriod(scalex)=mean_period;
    forcedBurstDur(scalex)=mean_burst_dur;
    forcedSPB(scalex)=spikes_per_burst;

end

openRec=[openIBI' openPeriod' openBurstDur' openSPB'];
closedRec=[closedIBI' closedPeriod' closedBurstDur' closedSPB'];
forcedRec=[forcedIBI' forcedPeriod' forcedBurstDur' forcedSPB'];

dlmwrite('open_IBI_BurstDur_SPB.csv',openRec,'precision',10)
dlmwrite('closed_IBI_BurstDur_SPB.csv',closedRec,'precision',10)
dlmwrite('forced_IBI_BurstDur_SPB.csv',forcedRec,'precision',10)


%% make plot

figure(1)

lw=2;

subplot(3,1,1)
hold on
plot(scale_factors,openIBI,'bs-','Linewidth',lw,'MarkerFaceColor','b')
plot(scale_factors,closedIBI,'ko-','Linewidth',lw,'MarkerFaceColor','k')
plot(scale_factors,forcedIBI,'gs-','Linewidth',lw,'MarkerFaceColor','g')
ylabel('IBI (s)','Interpreter','Latex')
xlabel('Scale Factor $\gamma$','Interpreter','Latex')
xlim([0.8 4.5])
set(gca,'box','off','XTick',[1:4],'YTick',[0:5:30])
h=legend(' $\bar{\tau}_h = \gamma \times 10^4 \quad \mbox{(open loop)}$',' $\bar{\tau}_h = \gamma \times 10^4 \quad \mbox{(closed loop)}$',' $\tau_{P_\mathrm{a}\mathrm{O}_2} =  \gamma$','Location','Northwest');
legend('boxoff')
set(h,'Interpreter','latex','fontsize',24)
grid on

subplot(3,1,2)
hold on
plot(scale_factors,openBurstDur,'bs-','Linewidth',lw,'MarkerFaceColor','b')
plot(scale_factors,closedBurstDur,'ko-','Linewidth',lw,'MarkerFaceColor','k')
plot(scale_factors,forcedBurstDur,'gs-','Linewidth',lw,'MarkerFaceColor','g')
ylabel('Burst Duration (s)','Interpreter','Latex')
xlabel('Scale Factor $\gamma$','Interpreter','Latex')
xlim([0.8 4.5])
set(gca,'box','off','XTick',[1:4],'YTick',[0:.5:2.5])
grid on

subplot(3,1,3)
hold on
plot(scale_factors,openSPB,'bs-','Linewidth',lw,'MarkerFaceColor','b')
plot(scale_factors,closedSPB,'ko-','Linewidth',lw,'MarkerFaceColor','k')
plot(scale_factors,forcedSPB,'gs-','Linewidth',lw,'MarkerFaceColor','g')
ylabel('Spikes Per Burst','Interpreter','Latex')
xlabel('Scale Factor $\gamma$','Interpreter','Latex')
xlim([0.8 4.5])
set(gca,'box','off','XTick',[1:4],'YTick',[0:10:60])
ylim([0 60])
grid on




