
function [ vsim, t, spikes_finalized, m_out, n_out, h_out, p_out ] = fh_sim2(dt, tmax, Iapp)
% Frankenhaeuser-Huxley model 
% returns simulated membrane potential vector (vsim) 
% and spike train vector (spikes_finalized)
% spike train vector indicates "time points" where spikes occured
% Iapp is pulse input current vector, generated from pulsemaker.m
% Iapp time should be equal to tmax

function expM1 = computeExpM1(x, y)
    if abs(x/y) < 1e-6
        expM1 = y*(1.0 - x/y/2.0);
    else
        expM1 = x/(exp(x/y)-1.0);
    end
end

function ghk = computeGhz(v, ci, co, tc, R, F, Vre)
    T = tc + 273.15;
    E = v + Vre;
    z = 1e-3 * F * E / (R * T);
    eco = co * efun(z);
    eci = ci * efun(-1.0*z);
    
    ghk = 0.001 * F * (eci - eco);
end

function efunn = efun(z)
    if abs(z) < 1.0e-4
        efunn = 1 - z/2.0;
    else
        efunn = z/(exp(z) - 1.0);
    end
end
    
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Constants

F=96484.6;              % Faraday's constant
R=8.3145;               % gas constant
celsius = 24;           % Temperature, in celsius

% Timestep and time vector 
% dt = 2e-8;              % Time step (s) must be extremely small
% tmax=0.35;              % Maximum simulation time
t=0:dt:tmax;            % Time vector

% Ion concentrations
C_o_Na=114.5;        % outside sodium concentration 
C_i_Na=13.74;        % inside sodium concentration
C_o_K=2.5;           % outside potassium concentration
C_i_K=120;           % inside potassium concentration

V_L=0.026;           % leak current potential mV (V_L)
Vre=-70;             % rest potential 
Cm=2e-6;             % Membrane capacitance (uF/cm^2)

% Permeabilities
na_per_con=8e-3;        % sodium permeability constant
k_per_con=0.54e-3;      % potassium permeability constant
nsp_per_con=1.2e-3;     % non-specific permeability constant

% Leak conductance 
G_L=0.0303; 

% FH Constants
% gate_m_alpha
A_alpha_m=0.36;    
B_alpha_m=22;      
C_alpha_m=3;       

% gate_m_beta
A_beta_m=0.4; 
B_beta_m=13; 
C_beta_m=20;

% gate_n_alpha
A_alpha_n=0.02; 
B_alpha_n=35; 
C_alpha_n=10;

% gate_n_beta
A_beta_n=0.05; 
B_beta_n=10; 
C_beta_n=10;

% gate_h_alpha
A_alpha_h=0.1;    % should be larger?
B_alpha_h=-10;   % should be smaller?
C_alpha_h=6;

% gate_h_beta
A_beta_h=4.5; 
B_beta_h=45; 
C_beta_h=10;

% gate_p_alpha
A_alpha_p=0.006; 
B_alpha_p=40; 
C_alpha_p=10;

% gate_p_beta
A_beta_p=0.09; 
B_beta_p=-25; 
C_beta_p=20;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initial Conditions

V0=0;              % Default initial condition for V
m0=0.0005;              % Default initial condition for m
h0=0.8249;              % Default initial condition for h
n0=0.0268;              % Default initial condition for n
p0=0.0049;              % Default initial condition for p

vsim=zeros(size(t));    % voltage vector
vsim(1)=V0;             % set the inititial value of voltage     
m=zeros(size(t));       % m: sodium activation gating variable
m(1)=m0;
n=zeros(size(t));       % n: potassium activation gating variable
n(1)=n0;  
h=zeros(size(t));       % h: sodim inactivation gating variable
h(1)=h0;               
p=zeros(size(t));       % p: non-specific permeability variable
p(1)=p0;

I_C=zeros(size(t));     % to plot and look at the total current
I_Na=zeros(size(t));    % to plot and look at sodium current
I_K=zeros(size(t));     % to plot and look at potassium current
I_p=zeros(size(t));     % to plot and look at non-specific delayed current
I_L=zeros(size(t));     % to plot and look at leak current

% Permeability vectors
P_Na=zeros(size(t));
P_K=zeros(size(t));
P_p=zeros(size(t));

% Spike detection
spikes = zeros(size(t)); % spike train vector
spike_flag = 0;          % for counting spikes. 0->1 if Vm > v_exceeds. 1->0 if Vm < v_unblock
v_exceeds = 70;       % threshold for detecting spikes
v_unblock = 10;      % threshold for after spike

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Simulation 

for i = 2:length(t) % now see how things change through time
    
    Vm = vsim(i-1);          % membrane potential for calculations
    Em = Vm + Vre;
    
    % First check if spikes occured
    if Vm > v_exceeds
        if spike_flag == 0
            spikes(i-1) = 1;     % we write 1 in spikes 
            spike_flag = 1;      % update spike flag
        end
    end
    
    if Vm < v_unblock
        spike_flag = 0;      % update spike flag
    end
    
    % Sodium and potassium gating variables are defined by the
    % voltage-dependent transition rates between states, labeled alpha and
    % beta.
    
    % sodium activation rate
    Q10Tpow = (celsius-20)*0.1;
    Q10T = 3.^Q10Tpow;
    
    alpha_m = A_alpha_m * computeExpM1((B_alpha_m-Vm), C_alpha_m);
    beta_m = A_beta_m * computeExpM1((Vm-B_beta_m), C_beta_m);
    % alpha_m = (A_alpha_m*(Vm-B_alpha_m))/(1-exp((B_alpha_m-Vm)/C_alpha_m));
    % beta_m = (A_beta_m*(B_beta_m-Vm))/(1-exp((Vm-B_beta_m)/C_beta_m));
    
    % potassium activation rate
    alpha_n = A_alpha_n * computeExpM1((B_alpha_n-Vm), C_alpha_n);
    beta_n = A_beta_n * computeExpM1((Vm-B_beta_n), C_beta_n);
    % alpha_n = A_alpha_n*(Vm-B_alpha_n)/(1-exp((B_alpha_n-Vm)/C_alpha_n));
    % beta_n = A_beta_n*(B_beta_n-Vm)/(1-exp((Vm-B_beta_n)/C_beta_n));
    
    % sodium inactivation rate 
    alpha_h = A_alpha_h * computeExpM1((Vm-B_alpha_h), C_alpha_h);
    % alpha_h = (A_alpha_h*(B_alpha_h-Vm))/(1-exp((Vm-B_alpha_h)/C_alpha_h));
    beta_h = A_beta_h/(1+exp((B_beta_h-Vm)/C_beta_h));
    
    % non-specific permeability rate
    alpha_p = A_alpha_p * computeExpM1((B_alpha_p-Vm), C_alpha_p);
    beta_p = A_beta_p * computeExpM1((Vm-B_beta_p), C_beta_p);
    % alpha_p = A_alpha_p*(Vm-B_alpha_p)/(1-exp((B_alpha_p-Vm)/C_alpha_p));
    % beta_p = A_beta_p*(B_beta_p-Vm)/(1-exp((Vm-B_beta_p)/C_beta_p));
   
    
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h,n & p.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    tau_p = 1/(alpha_p+beta_p);
    p_inf = alpha_p/(alpha_p+beta_p);
    
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt/tau_m;    % Update m
    
    h(i) = h(i-1) + (h_inf-h(i-1))*dt/tau_h;    % Update h
    
    n(i) = n(i-1) + (n_inf-n(i-1))*dt/tau_n;    % Update n
    
    p(i) = p(i-1) + (p_inf-p(i-1))*dt/tau_p;    % Update p
    
    % Calculate Permeability 
    P_Na(i) = na_per_con * h(i)*m(i)*m(i);        % Sodium permeability
    
    P_K(i) = k_per_con * n(i)*n(i);               % Potassium permeability
    
    P_p(i) = nsp_per_con * p(i)*p(i);             % nonspecific permeability
    
    
    % Update Current
    ZK = computeGhz(Vm, C_i_K, C_o_K, celsius, R, F, Vre);
    ZNa = computeGhz(Vm, C_i_Na, C_o_Na, celsius, R, F, Vre);
    
    I_Na(i) = P_Na(i) * ZNa;
    I_K(i) = P_K(i) * ZK;
    I_p(i) = P_p(i) * ZNa;
    I_L(i) = G_L * (Vm-V_L);
    I_C(i) = Iapp(i)-(I_Na(i)+I_K(i)+I_p(i)+I_L(i));
    vsim(i) = vsim(i-1) + I_C(i)*dt/Cm;
    
%     conc_na = C_o_Na - C_i_Na;
%     conc_k = C_o_K - C_i_K;
%     kelvin = celsius + 273.15;
%     efrt = Em*F / R*kelvin;
%     
%     I_Na(i) = P_Na(i) * F*efrt * (conc_na*exp(efrt)/(1-exp(efrt))); % total sodium current
%     
%     I_K(i) = P_K(i) * F*efrt * (conc_k*exp(efrt)/(1-exp(efrt))); % total potassium current
%     
%     I_p(i) = P_p(i) * F*efrt * (conc_na*exp(efrt)/(1-exp(efrt))); % total non-specific current
%     
%     I_L(i) = G_L*(Vm-V_L);    % Leak current is straightforward
%     
%     I_C(i) = Iapp(i)-(I_Na(i)+I_K(i)+I_p(i)+I_L(i)); % total current is sum of leak + active channels + applied current
%     
%     vsim(i) = vsim(i-1) + I_C(i)*dt/Cm;        % Update the membrane potential, V.
    
end

spikes_finalized = dt*find(spikes);             % Multiply dt to match actual time

m_out = m;                                      % Output m,n,h,p vector after simulation (10/12/22 update)
n_out = n;
h_out = h;
p_out = p;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
