
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Import Functions

import pulsemaker.* fh_sim2.*
import chop_train.* distance_matrix.* points.* ...
    information_from_matrix.* background.* golden_h.*

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Default plot settings

set(0,'DefaultLineLineWidth',1,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Part 1. Single spike experiment
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize essential inputs
dt=1e-7;
tmax=20;
n_trials=1;
window_len = 1;
tau = 0.3;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize and produce Iapp with pulsemaker

istart=2;
ibase=0;

npulses=5;
sep=2;
int=1;
len=0.12;

sigma_sep=0.2;
sigma_int=0;
sigma_len=0;

[ tvec, Iapp ] = pulsemaker(istart, ibase, npulses, sep, int, len, ...
    sigma_sep, sigma_int, sigma_len, tmax, dt);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Simulation 
[ V, t, spikes, m, n, h, p ] = fh_sim2(dt, tmax, Iapp);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot results - with m,n,h,p change

% input current
subplot(2,1,1);
plot(tvec, Iapp, 'g');
xlabel('Time (ms)')
ylabel('Input current (mA/cm^2)')
ylim([-1, 3]);
xlim([0 16]);

hold on;
subplot(2,1,2);
plot(t, V, 'g');
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)');
ylim([-10 120]);
xlim([0 16]);

% hold on;
% subplot(6,1,3);
% plot(t, m, 'r');
% xlabel('Time (s)')
% ylabel('m');
% 
% hold on;
% subplot(6,1,4);
% plot(t, n, 'g');
% xlabel('Time (s)')
% ylabel('n');
% 
% hold on;
% subplot(6,1,5);
% plot(t, h, 'b');
% xlabel('Time (s)')
% ylabel('h');
% 
% hold on;
% subplot(6,1,6);
% plot(t, p, 'm');
% xlabel('Time (s)')
% ylabel('p');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Part 2. Spike trains and MI calculation

% tmax = 20;
% istart = 0;     % Default time applied current starts
% ibase = 0;       % Default baseline current before/after pulse
% npulses = 1;    % Default number of current pulses
% 
% tau_sep = 2;     % mean of interval between pulses
% tau_int = 1;     % mean of pulse intensity 
% tau_len = 0.5;  % mean of one pulse length
% 
% % Simulate standard model
% sigma_sep = 0;
% sigma_int = 0;
% sigma_len = 0;
% 
% [ tvec, Iapp ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int, tau_len, ...
%     sigma_sep, sigma_int, sigma_len, tmax, dt);
% 
% % Simulation 
% [ V, t, spikes, mout, nout, hout, pout ] = fh_sim2(dt, tmax, Iapp, 24);
% 
% len_1 = length(spikes);
% frag_1 = chop_train(spikes, len_1, window_len, 0, tmax);   
% mat1 = distance_matrix(frag_1, tau);
% n1 = length(mat1);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #1: Input pulse interval & MI 
% Same intensity, same pulse length and same number of pulses
% Change the interval: choose the pulse interval from a normal distribution of 
% N(tau_sep, sigma_sep)
% Increase sigma_sep and observe MI change
% For each sigma_sep, average values of 'n_trials' trials

% sigma_sep_vec = 0:0.001:0.02;
% mi_vec = zeros(size(sigma_sep_vec));
% 
% for s=1:length(sigma_sep_vec)
%     info_vec = zeros(1, n_trials);
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = 10;
%         old_h = biggest_h;
% 
%         sigma_sepy = sigma_sep_vec(s);
% 
%         % Generate input currents
%         [ tvec, ivec ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, sigma_sepy, 0, 0, tmax, dt);
% 
%         % Simulate comparison model
%         [ fh_simulated, tt, spike_sim_fh, mm, nn, hh, pp ] = fh_sim2(dt, tmax, ivec);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_fh);
%         frag_2 = chop_train(spike_sim_fh, len_2, window_len, 0, tmax);   
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
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec(k) = correct_info;
%         
%     end
%     
%     mi_vec(s) = mean(info_vec);
%     
% end
% 
% % Plot results
% plot(sigma_sep_vec, mi_vec, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{sep}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #2: Input pulse intensity & MI 
% Same interval, same pulse length and same number of pulses
% Change the intensity: choose the pulse intensity from a normal distribution of 
% N(tau_int, sigma_int)
% Increase sigma_int and observe MI change
% For each sigma_int, average values of 'n_trials' trials

% sigma_int_vec = 0:0.01:0.2;
% mi_vec2 = zeros(size(sigma_int_vec));
% 
% for s=1:length(sigma_int_vec)
%     info_vec2 = zeros(1, n_trials);
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = 10;
%         old_h = biggest_h;
% 
%         sigma_inty = sigma_int_vec(s);
% 
%         % Generate input currents
%         [ tvec2, ivec2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, sigma_inty, 0, 0, tmax, dt);
% 
%         % Simulate comparison model
%         [ fh_simulated2, tt, spike_sim_fh2, mm, nn, hh, pp ] = fh_sim2(dt, tmax, ivec2);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_fh2);
%         frag_2 = chop_train(spike_sim_fh2, len_2, window_len, 0, tmax);   
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
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec2(k) = correct_info;
%         
%     end
%     
%     mi_vec2(s) = mean(info_vec2);
%     
% end
% 
% % Plot results
% plot(sigma_int_vec, mi_vec2, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{int}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Experiment #3: Input pulse length & MI 
% Same interval, same intensity and same number of pulses
% Change the pulse length: choose the pulse intensity from a normal distribution of 
% N(tau_len, sigma_len)
% Increase sigma_len and observe MI change
% For each sigma_len, average values of 'n_trials' trials

% sigma_len_vec = 0:0.01:0.2;
% mi_vec3 = zeros(size(sigma_len_vec));
% 
% for s=1:length(sigma_len_vec)
%     info_vec3 = zeros(1, n_trials);
%     for k=1:n_trials
%         % Upper bound on h
%         biggest_h = 10;
%         old_h = biggest_h;
% 
%         sigma_leny = sigma_len_vec(s);
% 
%         % Generate input currents
%         [ tvec2, ivec2 ] = pulsemaker(istart, ibase, npulses, tau_sep, tau_int,...
%         tau_len, sigma_leny, 0, 0, tmax, dt);
% 
%         % Simulate comparison model
%         [ fh_simulated3, tt, spike_sim_fh3, mm, nn, hh, pp ] = fh_sim2(dt, tmax, ivec2);
% 
%         % Chop spike trains
%         len_2 = length(spike_sim_fh3);
%         frag_2 = chop_train(spike_sim_fh3, len_2, window_len, 0, tmax);   
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
%         correct_info = MI - bias;
% 
%         % Save correct_info 
%         info_vec3(k) = correct_info;
%         
%     end
%     
%     mi_vec3(s) = mean(info_vec3);
%     
% end

% % Plot results
% plot(sigma_len_vec, mi_vec3, 'k');
% ylabel('MI (bits)');
% xlabel('sigma_{len}');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% subplot(2,1,1);
% plot(tvec, Iapp, 'k');
% xlabel('Time (ms)')
% ylabel('Input current (mA/cm^2)')
% ylim([-0.5, 1.5]);

% hold on;
% subplot(2,1,2);
% plot(t, V, 'k')
% xlabel('Time (ms)')
% ylabel('Membrane Potential (mV)');


