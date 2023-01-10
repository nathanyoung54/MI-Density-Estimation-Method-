% Experiment Codes for BC Senior Thesis
% Experiment #1: Plot F-I curve
% Experiment #2~#4: Validate MI estimator for pulse inputs

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Import Functions

% Experiment #1:
import hh_sim.* generate_iapp.* 

% Experiment #2~#4:
import pulsemaker.* chop_train.* distance_matrix.* points.* ...
    information_from_matrix.* background.* golden_h.*

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Default plot settings

set(0,'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define Parameters

global V_L E_Na E_K V0 G_L G_Na G_K Cm 
global tau_c S_bar_hh tau_r
global m0 h0 n0
global dt_hh dt_lif tmax

V_L = -0.070;       % Leak reversal potential (V)
E_Na = 0.045;       % Reversal for sodium channels (V)
E_K = -0.082;       % Reversal for potassium channels (V)
V0 = -0.065;        % Default initial condition for V

G_L = 30e-9;        % Leak conductance (S)
G_Na = 12e-6;       % Sodium conductance (S)
G_K = 3.6e-6;       % Potassium conductance (S)

Cm = 100e-12;       % Membrane capacitance (F)

tau_c = 10e-3;
S_bar_hh = 0.8e-9;  % S_bar = current. 30mV for voltage input
tau_r = 2e-3;       % refractory period

m0 = 0.05;          % Default initial condition for m
h0 = 0.5;           % Default initial condition for h
n0 = 0.35;          % Default initial condition for n

dt_hh = 1e-7;       % dt = 0.0001ms for H-H
dt_lif = 1e-4;      % dt = 0.1ms for LIF

tmax = 2;           % Maximum simulation time
tvec = 0:dt_hh:tmax;  % Time vector

window_len = 40e-3;     % time window for chopping trains
tau = 15e-3;
n_trials = 1;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize variables

vvec=zeros(size(tvec));    % voltage vector
vvec(1) = V0;              
nvec=zeros(size(tvec));    % n: potassium activation gating variable
nvec(1) = n0;              % initialize
mvec=zeros(size(tvec));    % m: sodium activation gating variable
mvec(1) = m0;              % initialize
hvec=zeros(size(tvec));    % h: sodim inactivation gating variable
hvec(1) = h0;              % initialize
        
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #1: F-I curve 

% ivalue = 0:1e-10:2e-9;                % make input current value vector
% fvec_hh = zeros(1,length(ivalue));       % initialize frequency vector
% 
% for i=1:length(ivalue)
%     curr_value = ivalue(i);
%     iapp = curr_value*ones(1,length(tvec));
%     
%     [ hh_simulated, spike_sim_hh ] = hh_sim(vvec, tvec, mvec, hvec, nvec, iapp);
%     
%     spike_count_hh = length(spike_sim_hh);
%     
%     freq_hh = spike_count_hh/tmax;
%     
%     fvec_hh(i) = freq_hh;
% end
% 
% % Plot Exp #1
% scatter(1e9*ivalue, fvec_hh, 'k', 'filled');
% ylabel('Firing rate (Hz)');
% xlabel('Applied Current (nA)');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #2: Input pulse interval & MI 
% Same intensity, same pulse length and same number of pulses
% Change the interval: choose the pulse interval from a normal distribution of 
% N(tau_sep, sigma_sep)
% Increase sigma_sep and observe MI change
% For each sigma_sep, average values of 'n_trials' trials

% Initiate variables and vectors (this holds for Experiment #3 and #4)
istart = 0.2;     % Default time applied current starts
ibase = 0e-9;        % Default baseline current before/after pulse
npulses = 40;        % Default number of current pulses

tau_sep = 20e-3;     % mean of interval between pulses
tau_int = 1e-9;    % mean of pulse intensity 
tau_len = 5e-3;      % mean of one pulse length

% Comment here to block Experiment #2

sigma_sep_vec = 0:1e-3:5e-3;
info_vec = zeros(1, n_trials);
mi_vec2 = zeros(size(sigma_sep_vec));


% Simulate standard model
sigma_sep = 0;
sigma_int = 0;
sigma_len = 0;

[ tvec2, ivec_stan ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
    tau_len, sigma_sep, sigma_int, sigma_len, tmax, dt_hh);


subplot(2,1,1);
plot(tvec2, 1e9*ivec_stan, 'b');
xlim([0 1.2]);
ylim([-1 3]);
ylabel('I_{app}(nA)')
xlabel('Time (s)');

[ hh_simulated, spike_sim_hh, mout, nout, hout ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec_stan);

subplot(2,1,2);
plot(tvec2, hh_simulated*1e3, 'b');
xlim([0 1.2]);
ylabel('Membrane Potential (mV)')
xlabel('Time (s)');
% 
% subplot(5,1,3);
% plot(tvec2, mout, 'k');
% xlim([0 1.2]);
% ylim([-0.2 1.2]);
% ylabel('m')
% xlabel('Time (s)');
% 
% subplot(5,1,4);
% plot(tvec2, nout, 'k');
% xlim([0 1.2]);
% ylim([-0.2 1.2]);
% ylabel('n')
% xlabel('Time (s)');
% 
% subplot(5,1,5);
% plot(tvec2, hout, 'k');
% xlim([0 1.2]);
% ylim([-0.2 1.2]);
% ylabel('h')
% xlabel('Time (s)');

% len_1 = length(spike_sim_hh);
% frag_1 = chop_train(spike_sim_hh, len_1, window_len, 0, tmax);   
% mat1 = distance_matrix(frag_1, tau);
% n1 = length(mat1);
% 
% for s=1:length(sigma_sep_vec)
%     
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = floor(2*tmax);
%         old_h = biggest_h;
% 
%         sigma_sepy = sigma_sep_vec(s);
% 
%         % Generate input currents
%         [ tvec2, ivec2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, sigma_sepy, 0, 0, tmax, dt_hh);
% 
%         % Simulate comparison model
%         [ hh_simulated2, spike_sim_hh2 ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec2);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_hh2);
%         frag_2 = chop_train(spike_sim_hh2, len_2, window_len, 0, tmax);   
% 
%         % Make distance matrices for each fragment
%         mat2 = distance_matrix(frag_2, tau);
%         n2 = length(mat2);
% 
%         % Golden Mean Search for h
%         [ h, old_h ] = golden_h(mat1, mat2, old_h, biggest_h);
% 
%         % Calculate MI with golden h
%         MI = information_from_matrix(mat1, mat2, h, h, 1);
%         bias = background(n1, h);
% 
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec(k) = correct_info;
%         
%     end
%     
%     mi_vec2(s) = mean(info_vec(k));
%     
% end
% 
% % Plot results
% plot(sigma_sep_vec, mi_vec2, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{sep}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #3: Input pulse intensity & MI 
% Same interval, same pulse length and same number of pulses
% Change the intensity: choose the pulse intensity from a normal distribution of 
% N(tau_int, sigma_int)
% Increase sigma_int and observe MI change
% For each sigma_int, average values of 'n_trials' trials

% Comment here to block Experiment #3
% 
% sigma_int_vec = 0:2e-10:1e-9;
% info_vec2 = zeros(1, n_trials);
% mi_vec3 = zeros(size(sigma_int_vec));
% 
% % Simulate standard model
% sigma_int = 0e-9;
% 
% [ tvec2, ivec_stan2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%     tau_len, 0, sigma_int, 0, tmax, dt_hh);
% 
% [ hh_simulated, spike_sim_hh ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec_stan2);
% 
% len_1 = length(spike_sim_hh);
% frag_1 = chop_train(spike_sim_hh, len_1, window_len, 0, tmax);   
% mat1 = distance_matrix(frag_1, tau);
% n1 = length(mat1);
% 
% for s=1:length(sigma_int_vec)
%     
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = floor(2*tmax);
%         old_h = biggest_h;
% 
%         sigma_inty = sigma_int_vec(s);
% 
%         % Generate input currents
%         [ tvec2, ivec2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, 0, sigma_inty, 0, tmax, dt_hh);
% 
%         % Simulate comparison model
%         [ hh_simulated3, spike_sim_hh3 ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec2);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_hh3);
%         frag_2 = chop_train(spike_sim_hh3, len_2, window_len, 0, tmax);   
% 
%         % Make distance matrices for each fragment
%         mat2 = distance_matrix(frag_2, tau);
%         n2 = length(mat2);
% 
%         % Golden Mean Search for h
%         [ h, old_h ] = golden_h(mat1, mat2, old_h, biggest_h);
% 
%         % Calculate MI with golden h
%         MI = information_from_matrix(mat1, mat2, h, h, 1);
%         bias = background(n1, h);
% 
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec2(k) = correct_info;
%         
%     end
%     
%     mi_vec3(s) = mean(info_vec2(k));
%     
% end
% 
% % Plot results
% plot(sigma_int_vec, mi_vec3, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{int}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #4: Input pulse length & MI 
% Same interval, same intensity and same number of pulses
% Change the pulse length: choose the pulse intensity from a normal distribution of 
% N(tau_len, sigma_len)
% Increase sigma_len and observe MI change
% For each sigma_len, average values of 'n_trials' trials
% 
% Comment here to block Experiment #4

% sigma_len_vec = 0:2e-3:10e-3;
% info_vec3 = zeros(1, n_trials);
% mi_vec4 = zeros(size(sigma_len_vec));
% 
% % Simulate standard model
% sigma_len = 0;
% 
% [ tvec2, ivec_stan3 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%     tau_len, 0, 0, sigma_len, tmax, dt_hh);
% 
% [ hh_simulated, spike_sim_hh ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec_stan3);
% 
% len_1 = length(spike_sim_hh);
% frag_1 = chop_train(spike_sim_hh, len_1, window_len, 0, tmax);   
% mat1 = distance_matrix(frag_1, tau);
% n1 = length(mat1);
% 
% for s=1:length(sigma_len_vec)
%     
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = floor(2*tmax);
%         old_h = biggest_h;
% 
%         sigma_leny = sigma_len_vec(s);
% 
%         % Generate input currents
%         [ tvec2, ivec2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, 0, 0, sigma_leny, tmax, dt_hh);
% 
%         % Simulate comparison model
%         [ hh_simulated4, spike_sim_hh4 ] = hh_sim(vvec, tvec2, mvec, hvec, nvec, ivec2);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_hh4);
%         frag_2 = chop_train(spike_sim_hh4, len_2, window_len, 0, tmax);   
% 
%         % Make distance matrices for each fragment
%         mat2 = distance_matrix(frag_2, tau);
%         n2 = length(mat2);
% 
%         % Golden Mean Search for h
%         [ h, old_h ] = golden_h(mat1, mat2, old_h, biggest_h);
% 
%         % Calculate MI with golden h
%         MI = information_from_matrix(mat1, mat2, h, h, 1);
%         bias = background(n1, h);
% 
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec3(k) = correct_info;
%         
%     end
%     
%     mi_vec4(s) = mean(info_vec3(k));
%     
% end
% 
% % Plot results
% plot(sigma_len_vec, mi_vec4, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{len}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


