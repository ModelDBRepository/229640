% This script reads in pre-computed data files and recreates Figure 14. To
% run the simulations that create the pre-computed data files, see
% "runsims_figure14.m"

clear all

breakVals=0:.01:.6;
breakDurs=1000:1000:60000;

for ix=1:length(breakVals)
    for jx=1:length(breakDurs)
        
        breakVal=breakVals(ix);
        breakDur=breakDurs(jx);
        
        avgpo2blood(ix,jx)=dlmread(sprintf('data_figure14/avgpo2blood_breakVal%1.4f_breakDur%5.4f_12_29_15.csv',breakVal,breakDur));

    end
end


%% Make plot

set(0,'DefaultAxesFontSize',16)

figure(1)
pcolor(avgpo2blood)
colormap(parula)
h=colorbar  
set(gca,'XTick',5:10:length(breakDurs),'XTickLabel',breakDurs(5:10:end)/1000,'YTick',1:10:length(breakVals),'YTickLabel',breakVals(1:10:end))
%xlabel('Clamp duration (s)','Interpreter','latex')
xlabel('Duration of chemosensory interruption (s)','Interpreter','latex')
%ylabel('Clamp severity ($g_\mathrm{tonic}$)','Interpreter','latex')
ylabel('Value of $g_\mathrm{tonic}$ during interruption','Interpreter','latex')
ylabel(h,'$P_\mathrm{a}\mathrm{O}_2$ after chemosensory restoration','Interpreter','latex','FontSize',20)


