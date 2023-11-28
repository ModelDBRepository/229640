function z = fastsubsystem(t,u)

global po2blood taulb R Temp

v=u(1);
n=u(2);
h=u(3);
alpha=u(4);
vollung=u(5);
po2lung=u(6);

gnap=2.8; gna=28; gk=11.2; gl=2.8; 
Ena=50; Ek=-85; El=-65; Esyn=0;

r=0.001; Tmax=1; VT=2; Kp=5;
E1=0.0025; E2=0.4;

PO2ext=(760-47)*.21; Vol0=2;

gtonic=0.3*(1-tanh((po2blood-85)/30));

Jlb=(1/taulb)*(po2lung-po2blood)*(vollung/(R*Temp));

%% Parameter Values (from Best et al, 2005)

C=21;  %pF

theta_mp=-40;  %mV
sigma_mp=-6;  %mV
taumax_h=10000;  %ms
theta_h=-48;  %mV
sigma_h=6;  %mV

theta_m=-34;  %mV
sigma_m=-5;  %mV

taumax_n=10;  %ms
theta_n=-29; %mV
sigma_n=-4;  %mV


%% Steady state values and time constants
mp_inf=1/(1+exp((v-theta_mp)/sigma_mp));

m_inf=1/(1+exp((v-theta_m)/sigma_m));

h_inf=1/(1+exp((v-theta_h)/sigma_h));
tau_h=taumax_h/cosh((v-theta_h)/(2*sigma_h));

n_inf=1/(1+exp((v-theta_n)/sigma_n));
tau_n=taumax_n/cosh((v-theta_n)/(2*sigma_n));


%% Currents

Inap=gnap*mp_inf*h*(v-Ena);
Ina=gna*(m_inf^3)*(1-n)*(v-Ena);
Ik=gk*(n^4)*(v-Ek);
Il=gl*(v-El);
Itonic=gtonic*(v-Esyn);

NT=Tmax/(1+exp(-(v-VT)/Kp));

%% Lung Volume

dvolrhs=max(0,-E1*(vollung-Vol0)+E2*alpha);


%% Equations

z(1)=(-Inap-Ina-Ik-Il-Itonic)/C;
z(2)=(n_inf-n)/tau_n;
z(3)=(h_inf-h)/tau_h;
z(4)=r*NT*(1-alpha)-r*alpha;
z(5)=-E1*(vollung-Vol0)+E2*alpha;
z(6)=(1/vollung)*(PO2ext-po2lung)*dvolrhs-Jlb*(R*Temp/vollung);

z=z';








