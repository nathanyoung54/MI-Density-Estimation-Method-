% Main code for calculating MI between two LIF neurons through
% density estimation method
% Modified code based on Houghton 2019
% F-I curve of LIF model added (5/25/22)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Import Functions

import lif_sim.* generate_iapp.* chop_train.* distance_matrix.* ...
    points.* information_from_matrix.* background.*

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define Parameters

global EL Rm Cm Vth Vre EK S_bar tau_m tau_c tau_r 

EL = -70e-3;                                % leak potential
Rm = 100e6;                                 % resistance
Cm = 0.1e-9;                                % membrane capacitance
Vth = -60e-3;                               % threshold potential
Vre = -65e-3;                               % reset potential
EK = -80e-3;                                % Potassium Nernst potential
S_bar = 30e-3;

tau_m = 12e-3;                              % time constant, Cm/GL
tau_c = 30e-3;
tau_r = 2e-3;                               % refractory period

global dt_lif tmax

dt_lif = 0.0001;                            % dt = 0.1ms
tmax = 2;                                   

% Variables

started = 0;
ended = 1;                                % total time = ended - started
mu = 1;                                   % decides dependency of inputs
tau = 15e-3;
window_len = 45e-3;                         % length of spike train fragment

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Default plot settings

set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% F-I curve 

% tvecy = 0:dt_lif:tmax;                      % time vector
% ivalue = 0:1e-4:50e-3;                      % make input current value vector
% fvec_lif = zeros(1,length(ivalue));         % initialize frequency vector
% 
% for i=1:length(ivalue)
%     curr_value = ivalue(i);
%     iapp = curr_value*ones(1,length(tvecy));
%     
%     [ lif_simulated, spike_sim_lif ] = lif_sim(tvecy, tau_r, iapp);
%     
%     spike_count_lif = length(spike_sim_lif);
%     
%     freq_lif = spike_count_lif/tmax;
%     
%     fvec_lif(i) = freq_lif;
% end
% 
% % Plot F-I curve
% plot(1e9*ivalue/Rm, fvec_lif, 'k');
% ylabel('Firing rate (Hz)');
% xlabel('Applied Current (nA)');
%title('LIF');
%ylim([-5 90]);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate input currents

[ i1vec, i2vec, tvec ] = generate_iapp(started, ended, dt_lif, mu, S_bar, tau_c);

% Simulate LIF model

[ lif_sim1, spike_sim1 ] = lif_sim(tvec, tau_r, i1vec);
[ lif_sim2, spike_sim2 ] = lif_sim(tvec, tau_r, i2vec);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% To see actual trains 

% subplot(3,1,1);
% plot(tvec, 1e3*p1vec, 'b');
% hold on;
% plot(tvec, 1e3*s1vec, 'r');
% %xlabel('Time (s)')
% ylabel('Input voltage (mV)')
% legend('P', 'S');
% 
% subplot(3,1,2);
% plot(tvec, 1e3*i1vec, 'k');
% %xlabel('Time (s)')
% ylabel('Input voltage (mV)')
% 
% subplot(3,1,3)
% plot(tvec, 1e3*lif_sim1, 'k')
% %title('Train #1')
% ylabel('Membrane potential (mV)');
% xlabel('Time (s)');

% For visualizing lif output
plot(tvec, 1e3*lif_sim2, 'k')
ylabel('Membrane Potential (mV)');
xlabel('Time (s)');

subplot(2,2,1);
plot(tvec, 1e3*i1vec, 'k')
xlabel('Time (s)')
ylabel('Input voltage (mV)')

subplot(2,2,2);
plot(tvec, 1e3*i2vec, 'k')
xlabel('Time (s)')
ylabel('Input voltage (mV)')

subplot(2,2,3);
plot(tvec, 1e3*lif_sim1, 'k')
xlabel('Time (s)')
ylabel('Membrane Potential (mV)')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Check: Chop spike trains

len_1 = length(spike_sim1);
len_2 = length(spike_sim2);

frag_1 = chop_train(spike_sim1, len_1, window_len, started, ended);     % fragments for spike train 1 (u)
frag_2 = chop_train(spike_sim2, len_2, window_len, started, ended);     % fragments for spike train 2 (v)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Check: Make distance matrices for each fragment

mat1 = distance_matrix(frag_1, tau);
mat2 = distance_matrix(frag_2, tau);

n1 = length(mat1);
n2 = length(mat2);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Golden Mean Search for h

biggest_h = floor(2*ended);
old_h = biggest_h;
phi = (1.0 + sqrt(5.0))/2.0;
stride = min([2*old_h, n1, biggest_h]);

a = 10;
b = stride;

c = floor(b-(b-a)/phi);
d = floor(a+(b-a)/phi);

info_c = information_from_matrix(mat1, mat2, c, c, 1) - background(n1, c);
info_d = information_from_matrix(mat1, mat2, d, d, 1) - background(n1, d);

while abs(d-c)>2
    
    if info_c > info_d
        b = d;
        d = c;
        c = floor(b-(b-a)/phi);
        info_d = info_c;
        info_c = information_from_matrix(mat1, mat2, c, c, 1) - background(n1, c);
    else
        a = c;
        c = d;
        d = floor(a+(b-a)/phi);
        info_c = info_d;
        info_d = information_from_matrix(mat1, mat2, d, d, 1) - background(n1, d);

    end
end

h = floor((a+b)/2);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Check: Calculate MI

MI = information_from_matrix(mat1, mat2, h, h, 1);
bias = background(n1, h);

correct_info = MI - bias;

formatSpec = 'The Mutual Information between Train #1 and Train #2 is %2.5f bits \n';
fprintf(formatSpec, correct_info);
fprintf("Calculated by the density estimation method \n");

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++









