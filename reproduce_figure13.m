clear all

%% Simulate trajectories for Panel A

initsR=[-59.166808240000002 5.302126437000001e-04 0.865394581200000 4.865298410000000e-06 2.000778444000000 42.071422380000001 41.832960610000001];
initsF=[-59.166808230000001 5.302126443000000e-04 0.865394581700000 4.865298416000000e-06 2.000778444000000 42.071374689999999 41.832913580000003];

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@events_figure13);

t0=0;
tf=6e4;

[tR,uR,teR,yeR] = ode15s(@closedloop,[0 tf],initsR,options);
[tF,uF,teF,yeF] = ode15s(@closedloop,[0 tf],initsF,options);

hR=uR(:,3);
volR=uR(:,5);
PO2bloodR=uR(:,7);

hF=uF(:,3);
volF=uF(:,5);
PO2bloodF=uF(:,7);

period=teR(7)-teR(6); % one period of boundary limit cycle

figure(100)
plot(tR,PO2bloodR)
xlabel('tR')
ylabel('PO2bloodR')
title('fig14-1-31-17')

figure(101)
plot(tF,PO2bloodF)
xlabel('tF')
ylabel('PO2bloodF')
title('fig14-1-31-17')

figure(102)
hold on
plot(tR,volR,'linewidth',2)
plot(teR(6),yeR(6,5),'ro')

figure(103)
hold on
plot(tR,uR(:,4),'g','linewidth',2)
plot(teR(6),yeR(6,4),'ro')

%% point b

figure(104)
solby=solb.y;
volb=solby(5,:);
plot(solb.x,volb,'b','linewidth',2)

figure(105)
solby=solb.y;
alphab=solby(4,:);
plot(solb.x,alphab,'g','linewidth',2)

%% point c

figure(106)
solcy=solc.y;
volc=solcy(5,:);
plot(solc.x,volc,'b','linewidth',2)

figure(107)
solcy=solc.y;
alphac=solcy(4,:);
plot(solc.x,alphac,'g','linewidth',2)

%% point d

figure(108)
soldy=sold.y;
vold=soldy(5,:);
plot(sold.x,vold,'b','linewidth',2)

figure(109)
soldy=sold.y;
alphad=soldy(4,:);
plot(sold.x,alphad,'g','linewidth',2)

%% point e

figure(110)
soley=sole.y;
vole=soley(5,:);
plot(sole.x,vole,'b','linewidth',2)

figure(111)
soley=sole.y;
alphae=soley(4,:);
plot(sole.x,alphae,'g','linewidth',2)

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Floquet eigenvectors for arrows and panels B-E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% arrow b and Panel B

initsb=yeR(6,:); % basepoint in quiescent phase

% simulate one unperturbed period
solb=ode15s(@closedloop,[0 period],initsb,options);

% make perturbations and simulate one period
eps=1e-7;

vb=initsb(1); nb=initsb(2); hb=initsb(3); alphab=initsb(4); vollungb=initsb(5); po2lungb=initsb(6); po2bloodb=initsb(7);

sol1b=ode15s(@closedloop,[0 period],[vb+eps nb hb alphab vollungb po2lungb po2bloodb],options);
sol2b=ode15s(@closedloop,[0 period],[vb nb+eps hb alphab vollungb po2lungb po2bloodb],options);
sol3b=ode15s(@closedloop,[0 period],[vb nb hb+eps alphab vollungb po2lungb po2bloodb],options);
sol4b=ode15s(@closedloop,[0 period],[vb nb hb alphab+eps vollungb po2lungb po2bloodb],options);
sol5b=ode15s(@closedloop,[0 period],[vb nb hb alphab vollungb+eps po2lungb po2bloodb],options);
sol6b=ode15s(@closedloop,[0 period],[vb nb hb alphab vollungb po2lungb+eps po2bloodb],options);
sol7b=ode15s(@closedloop,[0 period],[vb nb hb alphab vollungb po2lungb po2bloodb+eps],options);

% columns of multiplier matrix
x1b=sol1b.y(:,end)-solb.y(:,end);
x2b=sol2b.y(:,end)-solb.y(:,end);
x3b=sol3b.y(:,end)-solb.y(:,end);
x4b=sol4b.y(:,end)-solb.y(:,end);
x5b=sol5b.y(:,end)-solb.y(:,end);
x6b=sol6b.y(:,end)-solb.y(:,end);
x7b=sol7b.y(:,end)-solb.y(:,end);

Eb=[x1b x2b x3b x4b x5b x6b x7b]; % multiplier matrix

[evecB,evalB]=eig(Eb/eps); % eigenvalues of multiplier matrix are the Floquet multipliers

FloquetMultipliersB=diag(evalB)

% rescale eigenvectors for panel B

varMax=max(solb.y'); varMin=min(solb.y'); varDiff=abs(varMax-varMin);

s=varDiff; S=diag(s);

zeta1=inv(S)*evecB(:,1); zeta2=inv(S)*evecB(:,2); zeta3=inv(S)*evecB(:,3); zeta4=inv(S)*evecB(:,4); zeta5=inv(S)*evecB(:,5); zeta6=inv(S)*evecB(:,6); zeta7=inv(S)*evecB(:,7);

zeta1norm=zeta1/norm(zeta1); zeta2norm=zeta2/norm(zeta2); zeta3norm=zeta3/norm(zeta3); zeta4norm=zeta4/norm(zeta4); zeta5norm=zeta5/norm(zeta5); zeta6norm=zeta6/norm(zeta6); zeta7norm=zeta7/norm(zeta7);

zetaB=[zeta1norm zeta2norm zeta3norm zeta4norm zeta5norm zeta6norm zeta7norm];

%% arrow c and Panel C

% advance a fraction of one period to get basepoint shortly before first spike
[t1,u1,te1,ye1]=ode15s(@closedloop,[0 period*.2],initsb,options);

initsc=u1(end,:); % location c

% simulate one unperturbed period
solc=ode15s(@closedloop,[0 period],initsc,options);

% make perturbations and simulate one period
eps=1e-7;

vc=initsc(1); nc=initsc(2); hc=initsc(3); alphac=initsc(4); vollungc=initsc(5); po2lungc=initsc(6); po2bloodc=initsc(7);

sol1c=ode15s(@closedloop,[0 period],[vc+eps nc hc alphac vollungc po2lungc po2bloodc],options);
sol2c=ode15s(@closedloop,[0 period],[vc nc+eps hc alphac vollungc po2lungc po2bloodc],options);
sol3c=ode15s(@closedloop,[0 period],[vc nc hc+eps alphac vollungc po2lungc po2bloodc],options);
sol4c=ode15s(@closedloop,[0 period],[vc nc hc alphac+eps vollungc po2lungc po2bloodc],options);
sol5c=ode15s(@closedloop,[0 period],[vc nc hc alphac vollungc+eps po2lungc po2bloodc],options);
sol6c=ode15s(@closedloop,[0 period],[vc nc hc alphac vollungc po2lungc+eps po2bloodc],options);
sol7c=ode15s(@closedloop,[0 period],[vc nc hc alphac vollungc po2lungc po2bloodc+eps],options);

% columns of multiplier matrix
x1c=sol1c.y(:,end)-solc.y(:,end);
x2c=sol2c.y(:,end)-solc.y(:,end);
x3c=sol3c.y(:,end)-solc.y(:,end);
x4c=sol4c.y(:,end)-solc.y(:,end);
x5c=sol5c.y(:,end)-solc.y(:,end);
x6c=sol6c.y(:,end)-solc.y(:,end);
x7c=sol7c.y(:,end)-solc.y(:,end);

Ec=[x1c x2c x3c x4c x5c x6c x7c]; % multiplier matrix

[evecC,evalC]=eig(Ec/eps); % eigenvalues of multiplier matrix are the Floquet multipliers

FloquetMultipliersC=diag(evalC)

% rescale eigenvectors for panel C

varMax=max(solc.y'); varMin=min(solc.y'); varDiff=abs(varMax-varMin);

s=varDiff; S=diag(s);

zeta1=inv(S)*evecC(:,1); zeta2=inv(S)*evecC(:,2); zeta3=inv(S)*evecC(:,3); zeta4=inv(S)*evecC(:,4); zeta5=inv(S)*evecC(:,5); zeta6=inv(S)*evecC(:,6); zeta7=inv(S)*evecC(:,7);

zeta1norm=zeta1/norm(zeta1); zeta2norm=zeta2/norm(zeta2); zeta3norm=zeta3/norm(zeta3); zeta4norm=zeta4/norm(zeta4); zeta5norm=zeta5/norm(zeta5); zeta6norm=zeta6/norm(zeta6); zeta7norm=zeta7/norm(zeta7);

zetaC=[zeta1norm zeta2norm zeta3norm zeta4norm zeta5norm zeta6norm zeta7norm];

%% arrow d and Panel D

% advance a fraction of one period to get basepoint during active phase
[t1,u1,te1,ye1]=ode15s(@closedloop,[0 period*.35],initsb,options);

initsd=u1(end,:); % location d

% simulate one unperturbed period
sold=ode15s(@closedloop,[0 period],initsd,options);

% make perturbations and simulate one period
eps=1e-7;

vd=initsd(1); nd=initsd(2); hd=initsd(3); alphad=initsd(4); vollungd=initsd(5); po2lungd=initsd(6); po2bloodd=initsd(7);

sol1d=ode15s(@closedloop,[0 period],[vd+eps nd hd alphad vollungd po2lungd po2bloodd],options);
sol2d=ode15s(@closedloop,[0 period],[vd nd+eps hd alphad vollungd po2lungd po2bloodd],options);
sol3d=ode15s(@closedloop,[0 period],[vd nd hd+eps alphad vollungd po2lungd po2bloodd],options);
sol4d=ode15s(@closedloop,[0 period],[vd nd hd alphad+eps vollungd po2lungd po2bloodd],options);
sol5d=ode15s(@closedloop,[0 period],[vd nd hd alphad vollungd+eps po2lungd po2bloodd],options);
sol6d=ode15s(@closedloop,[0 period],[vd nd hd alphad vollungd po2lungd+eps po2bloodd],options);
sol7d=ode15s(@closedloop,[0 period],[vd nd hd alphad vollungd po2lungd po2bloodd+eps],options);

% columns of multiplier matrix
x1d=sol1d.y(:,end)-sold.y(:,end);
x2d=sol2d.y(:,end)-sold.y(:,end);
x3d=sol3d.y(:,end)-sold.y(:,end);
x4d=sol4d.y(:,end)-sold.y(:,end);
x5d=sol5d.y(:,end)-sold.y(:,end);
x6d=sol6d.y(:,end)-sold.y(:,end);
x7d=sol7d.y(:,end)-sold.y(:,end);

Ed=[x1d x2d x3d x4d x5d x6d x7d]; % multiplier matrix

[evecD,evalD]=eig(Ed/eps); % eigenvalues of multiplier matrix are the Floquet multipliers

FloquetMultipliersD=diag(evalD)

% rescale eigenvectors for panel D

varMax=max(sold.y'); varMin=min(sold.y'); varDiff=abs(varMax-varMin);

s=varDiff; S=diag(s);

zeta1=inv(S)*evecD(:,1); zeta2=inv(S)*evecD(:,2); zeta3=inv(S)*evecD(:,3); zeta4=inv(S)*evecD(:,4); zeta5=inv(S)*evecD(:,5); zeta6=inv(S)*evecD(:,6); zeta7=inv(S)*evecD(:,7);

zeta1norm=zeta1/norm(zeta1); zeta2norm=zeta2/norm(zeta2); zeta3norm=zeta3/norm(zeta3); zeta4norm=zeta4/norm(zeta4); zeta5norm=zeta5/norm(zeta5); zeta6norm=zeta6/norm(zeta6); zeta7norm=zeta7/norm(zeta7);

zetaD=[zeta1norm zeta2norm zeta3norm zeta4norm zeta5norm zeta6norm zeta7norm];

%% arrow e and Panel E

% advance a fraction of one period to get basepoint shortly after last spike
[t1,u1,te1,ye1]=ode15s(@closedloop,[0 period*.7],initsb,options);

initse=u1(end,:); % location e

% simulate one unperturbed period
sole=ode15s(@closedloop,[0 period],initse,options);

% make perturbations and simulate one period
eps=1e-7;

ve=initse(1); ne=initse(2); he=initse(3); alphae=initse(4); vollunge=initse(5); po2lunge=initse(6); po2bloode=initse(7);

sol1e=ode15s(@closedloop,[0 period],[ve+eps ne he alphae vollunge po2lunge po2bloode],options);
sol2e=ode15s(@closedloop,[0 period],[ve ne+eps he alphae vollunge po2lunge po2bloode],options);
sol3e=ode15s(@closedloop,[0 period],[ve ne he+eps alphae vollunge po2lunge po2bloode],options);
sol4e=ode15s(@closedloop,[0 period],[ve ne he alphae+eps vollunge po2lunge po2bloode],options);
sol5e=ode15s(@closedloop,[0 period],[ve ne he alphae vollunge+eps po2lunge po2bloode],options);
sol6e=ode15s(@closedloop,[0 period],[ve ne he alphae vollunge po2lunge+eps po2bloode],options);
sol7e=ode15s(@closedloop,[0 period],[ve ne he alphae vollunge po2lunge po2bloode+eps],options);

% columns of multiplier matrix
x1e=sol1e.y(:,end)-sole.y(:,end);
x2e=sol2e.y(:,end)-sole.y(:,end);
x3e=sol3e.y(:,end)-sole.y(:,end);
x4e=sol4e.y(:,end)-sole.y(:,end);
x5e=sol5e.y(:,end)-sole.y(:,end);
x6e=sol6e.y(:,end)-sole.y(:,end);
x7e=sol7e.y(:,end)-sole.y(:,end);

Ee=[x1e x2e x3e x4e x5e x6e x7e]; % multiplier matrix

[evecE,evalE]=eig(Ee/eps); % eigenvalues of multiplier matrix are the Floquet multipliers

FloquetMultipliersE=diag(evalE)

% rescale eigenvectors for panel E

varMax=max(sole.y'); varMin=min(sole.y'); varDiff=abs(varMax-varMin);

s=varDiff; S=diag(s);

zeta1=inv(S)*evecE(:,1); zeta2=inv(S)*evecE(:,2); zeta3=inv(S)*evecE(:,3); zeta4=inv(S)*evecE(:,4); zeta5=inv(S)*evecE(:,5); zeta6=inv(S)*evecE(:,6); zeta7=inv(S)*evecE(:,7);

zeta1norm=zeta1/norm(zeta1); zeta2norm=zeta2/norm(zeta2); zeta3norm=zeta3/norm(zeta3); zeta4norm=zeta4/norm(zeta4); zeta5norm=zeta5/norm(zeta5); zeta6norm=zeta6/norm(zeta6); zeta7norm=zeta7/norm(zeta7);

zetaE=[zeta1norm zeta2norm zeta3norm zeta4norm zeta5norm zeta6norm zeta7norm];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close('all')

%% panel A

set(0,'DefaultAxesFontSize',24)

xloR=3e4;
xhiR=5.5e4;
xloF=3e4;
xhiF=4.22e4;

indsR=find(tR>xloR & tR <xhiR);
indsF=find(tF>xloF & tF <xhiF);

figure(1)
hold on
plot3(hR(indsR),volR(indsR),PO2bloodR(indsR),'g','Linewidth',1)
plot3(hF(indsF),volF(indsF),PO2bloodF(indsF),'r','Linewidth',1)
xlabel('$h$','Interpreter','Latex')
ylabel('vol$_\mathrm{L}$','Interpreter','Latex')
zlabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
view(3)
az=-32; el=-6;
grid on
set(gca,'box','off','XTick',[.5 .51 .52],'XTickLabel',{'0.50' '0.51' '0.52'},'YTick',[2.15 2.25 2.35],'ZTick',[76 77 78],'ZTickLabel',{'76' '77' '78'})
xlim([.498 .52])
ylim([2.15 2.35])
zlim([75.8 78])


%% add arrows to panel A

arrowb=initsb([3 5 7])-0.55*[evecB(3,1) evecB(5,1) evecB(7,1)];
arrowc=initsc([3 5 7])+0.5*[evecC(3,1) evecC(5,1) evecC(7,1)];
arrowd=initsd([3 5 7])-1*[evecD(3,1) evecD(5,1) evecD(7,1)];
arrowe=initse([3 5 7])+0.75*[evecE(3,1) evecE(5,1) evecE(7,1)];

figure(1)
arrow(initsb([3 5 7]),arrowb,'Linewidth',5)
arrow(initsc([3 5 7]),arrowc,'Linewidth',5)
arrow(initsd([3 5 7]),arrowd,'Linewidth',5)
arrow(initse([3 5 7]),arrowe,'Linewidth',5)


%% panels B through E

set(0,'DefaultAxesFontSize',36)

% panel B
figure(2)
labels={'V','n','h','\alpha','vol_L','P_AO_2','P_aO_2'};
bar(zetaB(:,1))
title('quiescent phase','Interpreter','latex','Fontsize',36)
set(gca,'box','off','XTickLabel',labels,'YTick',-.1:.2:.7)
rotateXLabels(gca,90)
ylim([-.1 .75])

% panel C
figure(3)
labels={'$V$','$n$','$h$','$\alpha$','vol$_\mathrm{L}$','$P_\mathrm{A}\mathrm{O}_2$','$P_\mathrm{a}\mathrm{O}_2$'};
bar(-zetaC(:,1))
title('first spike','Interpreter','latex','Fontsize',36)
set(gca,'box','off','XTickLabels',labels,'TickLabelInterpreter','latex','YTick',-.1:.2:.7)
rotateXLabels(gca,90)
ylim([-.1 .75])

% panel D
figure(4)
labels={'$V$','$n$','$h$','$\alpha$','vol$_\mathrm{L}$','$P_\mathrm{A}\mathrm{O}_2$','$P_\mathrm{a}\mathrm{O}_2$'};
bar(zetaD(:,1))
title('active phase','Interpreter','latex','Fontsize',36)
set(gca,'box','off','XTickLabels',labels,'TickLabelInterpreter','latex','YTick',-.1:.2:.7)
rotateXLabels(gca,90)
ylim([-.1 .75])

% panel E
figure(5)
labels={'$V$','$n$','$h$','$\alpha$','vol$_\mathrm{L}$','$P_\mathrm{A}\mathrm{O}_2$','$P_\mathrm{a}\mathrm{O}_2$'};
bar(-zetaE(:,1))
title('last spike','Interpreter','latex','Fontsize',36)
set(gca,'box','off','XTickLabels',labels,'TickLabelInterpreter','latex','YTick',-.1:.2:.7)
rotateXLabels(gca,90)
ylim([-.1 .75])