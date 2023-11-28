% This script reproduces Figure 7 by simulating the averaged system for
% a range of M values and finding the fixed points (gbar=0)

clear all

%% Parameters

global po2blood taulb R Temp

R=62.364;
Temp=310;
taulb=500;

po2bloodvals=110:-1:10;

Mvals=2e-6:0.1e-6:18e-6;

%% Run fast subsystem (closed-loop with constant PO2 blood) over a range of PO2 blood values

for px=1:length(po2bloodvals)

    po2blood=po2bloodvals(px)
    
    po2bloodTimes10=round(10*po2blood);
    po2bloodTimes10plus1=round(10*(po2blood+.1));
    
    % Initial Conditions      
    v0=-55; n0=.001; h0=.75; alpha0=.0001; vollung0=2; po2lung0=po2blood;
    inits0=[v0 n0 h0 alpha0 vollung0 po2lung0];
        
    t0=0;
    tf=60000*5;

    options=[];
    [t0,u0]=ode15s('fastsubsystem',[t0 tf],inits0,options);
    
    inits=u0(end,:); 
    
    %% Determine whether CPG is quiescient, spiking, or bursting

    options=odeset('RelTol',1e-8,'AbsTol',1e-8);
    newSol=ode15s(@fastsubsystem,[0 30000],inits,options);
    
    t=newSol.x;
    
    v=newSol.y(1,:);
    vollung=newSol.y(5,:);
    po2lung=newSol.y(6,:);
    Jlb=(1/taulb)*(po2lung-po2blood).*(vollung/(R*Temp));
        
    numT=length(t);
    
    spikeInds=[];
    tlap=0;

    for ix=2:numT
        if v(ix)>0 && v(ix-1)<0
            if t(ix)>(tlap+1)
                spikeInds=[spikeInds ix];
                tlap=t(ix);
            end
        end
    end

    spikeTimes=t(spikeInds);
    numSpikes=length(spikeTimes);
    
    if numSpikes<=1
        numBursts=0;
    end
    
    if numSpikes>1
   
        ISI=spikeTimes(2:end)-spikeTimes(1:(end-1));
        
        largeISIinds=find(ISI>0.5*max(ISI));
        numBursts=length(largeISIinds);
    
        if numBursts>1
        
            largeISIs=ISI(largeISIinds);

            lastBurstStartInd=largeISIinds(end);
            secondToLastBurstStartInd=largeISIinds(end-1);

            lastBurstStartTime=spikeTimes(lastBurstStartInd+1);
            secondToLastBurstStartTime=spikeTimes(secondToLastBurstStartInd+1);

            lastBurstInd=find(t==lastBurstStartTime);
            secondToLastBurstInd=find(t==secondToLastBurstStartTime);
            
        end
    
    end
    
    %% Obtain one period of quiescent, spiking, or bursting solution
    
    if numBursts>1
        period=t(lastBurstInd)-t(secondToLastBurstInd);
        tSub=t(secondToLastBurstInd:lastBurstInd);
        JlbSub=Jlb(secondToLastBurstInd:lastBurstInd);
    end
    
    if numBursts<=1
        if numSpikes>1
            period=t(spikeInds(end))-t(spikeInds(end-1));
            tSub=t(spikeInds(end-1):spikeInds(end));
            JlbSub=Jlb(spikeInds(end-1):spikeInds(end));
        end
        
        if numSpikes<=1
            period=1;
            tSub=t(end);
            JlbSub=Jlb(end);
        end
    
    end
    
    %% Calculate gbar for a range of M values
    
    for mx=1:length(Mvals)

        M=Mvals(mx);
        Hb=150; volblood=5; 
        eta=Hb*1.36; gamma=volblood/22400; betaO2=0.03;

        c=2.5; K=26;
        SaO2=(po2blood^c)/(po2blood^c+K^c);
        partial=(c*po2blood^(c-1))*(1/(po2blood^c+K^c)-(po2blood^c)/((po2blood^c+K^c)^2));
        CaO2=eta*SaO2+betaO2*po2blood;
        Jbt=M*CaO2*gamma;

        dxdt=(JlbSub-Jbt)/(gamma*(betaO2+eta*partial));

        if numSpikes>1

            intFy=trapz(tSub,dxdt);
            gbar(px,mx)=intFy/period;

        end

        if numSpikes<=1
            gbar(px,mx)=dxdt;
        end
    
    end     
end

%% Find location of fixed points (gbar=0)

fprec=zeros(length(Mvals),3);

for mx=1:length(Mvals)

    fp=[];

    for px=2:length(po2bloodvals)

            if gbar(px-1,mx)>0
                if gbar(px,mx)<0
                    fp=[fp mean(po2bloodvals((px-1):px))];
                end
            end

            if gbar(px-1,mx)<0
                if gbar(px,mx)>0
                    fp=[fp mean(po2bloodvals((px-1):px))];
                end
            end


    end

    fprec(mx,1:length(fp))=fp';
    
end


%% Make plots

set(0,'DefaultAxesFontSize',24)

% figure 8A

figure(1)
hold on
plot(po2bloodvals,gbar(:,21),'g','Linewidth',3)
plot(po2bloodvals,gbar(:,61),'c','Linewidth',3)
plot(po2bloodvals,gbar(:,141),'m','Linewidth',3)
plot(po2bloodvals,0*po2bloodvals,'k--','LInewidth',3)
plot(fprec(21,1),0,'ko','MarkerSize',10,'MarkerFaceColor','g')
plot(fprec(21,2),0,'ko','MarkerSize',10,'MarkerFaceColor','g')
plot(fprec(21,3),0,'ko','MarkerSize',10,'MarkerFaceColor','g')
plot(fprec(61,1),0,'ko','MarkerSize',10,'MarkerFaceColor','c')
plot(fprec(61,2),0,'ko','MarkerSize',10,'MarkerFaceColor','c')
plot(fprec(61,3),0,'ko','MarkerSize',10,'MarkerFaceColor','c')
plot(fprec(141,1),0,'ko','MarkerSize',10,'MarkerFaceColor','m')
h=legend('$M = 0.4 \times 10^{-5}$','$M = 0.8 \times 10^{-5}$','$M = 1.6 \times 10^{-5}$','Location','Southwest');
legend('boxoff')
set(h,'Interpreter','latex','FontSize',20)
xlabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','latex','Fontsize',24)
ylabel('$\bar{g}$','interpreter','latex','FontSize',24)
ylim([-12.5e-3 3.5e-3])
xlim([10 100])
set(gca,'box','off','YTick',-12e-3:2e-3:4e-3)
grid on

% figure 8B

figure(2)
hold on
plot(Mvals,fprec(:,1),'o','Color',[1 .5 0],'MarkerFaceColor',[1 .5 0],'MarkerSize',6)
plot(Mvals,fprec(:,2),'o','Color',[1 .5 0]','MarkerFaceColor',[1 .5 0],'MarkerSize',6)
plot(Mvals,fprec(:,3),'o','Color',[1 .5 0],'MarkerFaceColor',[1 .5 0],'MarkerSize',6)
xlim([.2e-5 18e-6])
ylim([1 100])
h=legend('$\bar{g} = 0$','Location','Northeast');
set(h,'interpreter','latex','FontSize',24)
legend('boxoff')
xlabel('$M$','Interpreter','latex','FontSize',24)
ylabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','latex','Fontsize',24)
set(gca,'XTick',[0.4e-5:.4e-5:1.6e-5])
grid on












