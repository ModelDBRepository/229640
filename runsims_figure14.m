clear all

global M taulb
global breakFlag breakVal

mkdir('newdata_runsims_figure14')

%% Parameters

breakVals=0:.01:.6;

%% for gtonic=0.1
breakDurs=1000:1000:60000;

for ix=1:length(breakVals)
    for jx=1:length(breakDurs)
        
        breakVal=breakVals(ix)
        breakDur=breakDurs(jx)

taulb=500;
M=8e-6;

% Initial conditions [V n h alpha vollung PO2lung PO2blood]
inits2=[-58.625 6.0707e-04 0.71925 0.0012881 2.3288 1.0331e+02 1.0222e+02];

tf=6e4;

options=odeset('RelTol',1e-6,'AbsTol',1e-9);

breakFlag=1;

[t2,u2]=ode15s('closedloop_clampgtonic',[tf tf+breakDur],inits2,options);

inits3=u2(end,:);
breakFlag=0;

[t3,u3]=ode15s('closedloop_clampgtonic',[tf+breakDur tf+breakDur+2*tf],inits3,options);

inits4=u3(end,:);
[t4,u4]=ode15s('closedloop_clampgtonic',[tf+breakDur+2*tf tf+breakDur+3*tf],inits4,options);

inits5=u4(end,:);
[t5,u5]=ode15s('closedloop_clampgtonic',[tf+breakDur+3*tf:.1:tf+breakDur+3*tf+10000],inits5,options);

po2blood=u5(:,end);

avgpo2blood(ix,jx)=0.5*(min(po2blood)+max(po2blood));

dlmwrite(sprintf('newdata_runsims_figure14/avgpo2blood_breakVal%1.4f_breakDur%5.4f.csv',breakVal,breakDur),avgpo2blood(ix,jx),'precision',10)

    end
    
end


%% Make plot

set(0,'DefaultAxesFontSize',20)

figure(1)
pcolor(avgpo2blood)
colormap(parula)
colorbar  
set(gca,'XTick',5:10:length(breakDurs),'XTickLabel',breakDurs(5:10:end)/1000,'YTick',1:10:length(breakVals),'YTickLabel',breakVals(1:10:end))
xlabel('Clamp duration (s)','Interpreter','latex')
ylabel('Clamp severity ($g_\mathrm{tonic}$)','Interpreter','latex')





