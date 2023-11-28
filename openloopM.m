function z = openloopM(t,u)

global gtonic_open M

%% State variables

v = u(1); n = u(2); h = u(3); alpha = u(4); vollung = u(5); PO2lung = u(6); PO2blood = u(7);

%% CPG 

% capacitance
C = 21;  

% maximal conductances
gnap=2.8; gna=28; gk=11.2; gl=2.8;

% reversal potentials
Ena=50; Ek=-85; El=-65; Esyn=0;

% persistent sodium
theta_mp = -40;  sigma_mp = -6; 
theta_h = -48; sigma_h = 6; taumax_h = 10000;

mp_inf = 1/(1+exp((v-theta_mp)/sigma_mp));
h_inf = 1/(1+exp((v-theta_h)/sigma_h));
tau_h = taumax_h/cosh((v-theta_h)/(2*sigma_h));

Inap = gnap*mp_inf*h*(v-Ena);

% transient sodium
theta_m = -34; sigma_m = -5;

m_inf = 1/(1+exp((v-theta_m)/sigma_m));

Ina = gna*(m_inf^3)*(1-n)*(v-Ena);

% potassium
theta_n = -29; sigma_n = -4; taumax_n = 10;

Ik = gk*(n^4)*(v-Ek);

n_inf = 1/(1+exp((v-theta_n)/sigma_n));
tau_n = taumax_n/cosh((v-theta_n)/(2*sigma_n));

% leak
Il = gl*(v-El);

% tonic
Itonic = gtonic_open*(v-Esyn);

%% Motor pool

r = 0.001; Tmax = 1; VT = 2; Kp = 5;

NT = Tmax/(1+exp(-(v-VT)/Kp));

%% Lung volume

E1 = 0.0025; E2 = 0.4; Vol0 = 2;

dvolrhs=max(0,-E1*(vollung-Vol0)+E2*alpha);

%% Lung oxygen

PO2ext = (760-47)*.21;  R = 62.364; Temp = 310; 

taulb = 500;

%% Blood oxygen

Hb = 150; volblood = 5; eta = Hb*1.36; gamma = volblood/22400; betaO2 = 0.03;

c = 2.5; K = 26;
SaO2 = (PO2blood^c)/(PO2blood^c+K^c);
CaO2 = eta*SaO2+betaO2*PO2blood;
partial = (c*PO2blood^(c-1))*(1/(PO2blood^c+K^c)-(PO2blood^c)/((PO2blood^c+K^c)^2));

Jlb=(1/taulb)*(PO2lung-PO2blood)*(vollung/(R*Temp));
Jbt=M*CaO2*gamma;


%% Differential equations

z(1) = (-Inap-Ina-Ik-Il-Itonic)/C;
z(2) = (n_inf-n)/tau_n;
z(3) = (h_inf-h)/tau_h;
z(4) = r*NT*(1-alpha)-r*alpha;
z(5) = -E1*(vollung-Vol0)+E2*alpha;
z(6) = (1/vollung)*(PO2ext-PO2lung)*dvolrhs-Jlb*(R*Temp/vollung);
z(7) = (Jlb-Jbt)/(gamma*(betaO2+eta*partial));

z=z';




