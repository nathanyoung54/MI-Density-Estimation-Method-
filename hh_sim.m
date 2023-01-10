
function [ vsim, spikes_finalized, m_out, n_out, h_out ] = hh_sim(V, t, m, h, n, I)
% Hodgkin-Huxley model 
% returns simulated membrane potential vector and spike train vector
% spike train vector indicates "time points" where spikes occured
% input will be empty membrane potential vector (V), initialized with V_L
% t is time vector
% m, h are sodium activation, inactivation variables (gating variables)
% n is potassium activation variable (gating variable)
% I is input current vector, generated from generate_Iapp.m

global V_L E_Na E_K  G_L G_Na G_K Cm dt_hh

Itot=zeros(size(t));     % to plot and look at the total current
I_Na=zeros(size(t));     % to plot and look at sodium current
I_K=zeros(size(t));      % to plot and look at potassium current
I_L=zeros(size(t));      % to plot and look at leak current
spikes = zeros(size(t)); % spike train vector
    
vsim = zeros(size(V));   % membrane potential vector (return vector)
vsim(1) = V_L;           % initialize with leak potential

spike_flag = 0;          % for counting spikes. 0->1 if mp > v_exceeds. 1->0 if mp < v_unblock
v_exceeds = -10e-3;      % threshold for detecting spikes
v_unblock = -30e-3;      % threshold for after spike

for i = 2:length(t)      % now see how things change through time
    
     Vm = vsim(i-1);          % membrane potential for calculations
     
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
    
    % First, sodium activation rate
    if ( Vm == -0.045 )     % to avoid dividing zero by zero
        alpha_m = 1e3;      % value calculated analytically
    else
        alpha_m = (1e5*(-Vm-0.045))/(exp(100*(-Vm-0.045))-1);
    end
    beta_m = 4000*exp((-Vm-0.070)/0.018);   % Sodium deactivation rate
    alpha_h = 70*exp(50*(-Vm-0.070));       % Sodium inactivation rate
    beta_h = 1000/(1+exp(100*(-Vm-0.040))); % Sodium deinactivation rate
    
    % Second, potassium activation rate
    if ( Vm == -0.060)      % to avoid dividing by zero
        alpha_n = 100;      % value calculated analytically
    else                    % potassium activation rate
        alpha_n = (1e4*(-Vm-0.060))/(exp(100*(-Vm-0.060))-1);
    end
    beta_n = 125*exp((-Vm-0.070)/0.08);     % potassium deactivation rate
        
    % From the alpha and beta for each gating variable we find the steady
    % state values (_inf) and the time constants (tau_) for each m,h and n.
    
    tau_m = 1/(alpha_m+beta_m);
    m_inf = alpha_m/(alpha_m+beta_m);
    
    tau_h = 1/(alpha_h+beta_h);
    h_inf = alpha_h/(alpha_h+beta_h);
    
    tau_n = 1/(alpha_n+beta_n);
    n_inf = alpha_n/(alpha_n+beta_n);
    
    % Updates
    
    m(i) = m(i-1) + (m_inf-m(i-1))*dt_hh/tau_m;    % Update m
    
    h(i) = h(i-1) + (h_inf-h(i-1))*dt_hh/tau_h;    % Update h
    
    n(i) = n(i-1) + (n_inf-n(i-1))*dt_hh/tau_n;    % Update n
    
    I_Na(i) = G_Na*m(i)*m(i)*m(i)*h(i)*(E_Na-vsim(i-1)); % total sodium current
    
    I_K(i) = G_K*n(i)*n(i)*n(i)*n(i)*(E_K-vsim(i-1)); % total potassium current
    
    I_L(i) = G_L*(V_L-vsim(i-1));    % Leak current is straightforward
    
    Itot(i) = I(i)+(I_L(i)+I_Na(i)+I_K(i)); % total current is sum of leak + active channels + applied current
    
    vsim(i) = vsim(i-1) + Itot(i)*dt_hh/Cm;        % Update the membrane potential, vsim for output
    
end

spikes_finalized = dt_hh*find(spikes);             % Multiply dt to match actual time

m_out = m;                                         % Output m,n,h vector after simulation (10/12/22 update)
n_out = n;
h_out = h;

end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++