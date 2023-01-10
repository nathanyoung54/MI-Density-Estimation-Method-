% Main code for calculating MI between two Hodgkin Huxley neurons through
% Density Estimation Method
% Modified code based on Houghton 2019

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Import Functions

import hh_sim.* generate_iapp.* chop_train.* distance_matrix.* points.* ...
    information_from_matrix.* background.*

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define Parameters

global V_L E_Na E_K V0 G_L G_Na G_K Cm 
global tau_c S_bar_hh
global m0 h0 n0


V_L = -0.070;       % Leak reversal potential (V)
E_Na = 0.045;       % Reversal for sodium channels (V)
E_K = -0.080;       % Reversal for potassium channels (V)
V0 = -0.065;        % Default initial condition for V

G_L = 30e-9;        % Leak conductance (S)
G_Na = 12e-6;       % Sodium conductance (S)
G_K = 3.6e-6;       % Potassium conductance (S)

Cm = 100e-12;       % Membrane capacitance (F)

tau_c = 30e-3;
S_bar_hh = 0.8e-9;  % S_bar = current. 30mV for voltage input

m0 = 0.05;          % Default initial condition for m
h0 = 0.5;           % Default initial condition for h
n0 = 0.35;          % Default initial condition for n

global dt_hh

dt_hh = 1e-7;          % dt = 0.0000ms

% Variables

mu = 0.7;               % mu "variable"
started = 0;            % starting time point 
ended = 1;            % ending time point
window_len = 40e-3;     % time window for chopping trains
tau = 15e-3;
h = 10;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize variables

tvec = started:dt_hh:ended;

vvec=zeros(size(tvec));       % voltage vector
vvec(1) = V0;                 % set the inititial value of voltage

nvec=zeros(size(tvec));       % n: potassium activation gating variable
nvec(1) = n0;                 % initialize as zero
mvec=zeros(size(tvec));       % m: sodium activation gating variable
mvec(1) = m0;                 % initialize as zero
hvec=zeros(size(tvec));       % h: sodim inactivation gating variable
hvec(1) = h0;                 % initialize as zero

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Default plot settings

set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Check: Plot spike trains

% Make Iapp vectors. isim1 and isim2 are Iapp vectors for neuron u and v

[ isim1, isim2, tsim ] = generate_iapp(started, ended, dt_hh, mu, S_bar_hh, tau_c);

% subplot(2,2,1);
% plot(tsim, 1e9*isim1, 'k');
% xlabel('Time (s)');
% ylabel('Input Current (nA)');
% hold on;

subplot(2,2,2);
plot(tsim, 1e9*isim2, 'k');
xlabel('Time (s)');
ylabel('Input Current (nA)');
hold on;

% Run H-H simulation. vsim1 and vsim2 are mp vectors for neuron u and v st1
% and st2 vectors are spike train vectors for neuron u and v

[ vsim1, spike_sim1 ] = hh_sim(vvec, tsim, mvec, hvec, nvec, isim1);

% subplot(2,2,3);
% plot(tsim, 1e3*vsim1, 'k');
% xlabel('Time (s)');
% ylabel('Membrane Potental (mV)');
% hold on;


[ vsim2, spike_sim2 ] = hh_sim(vvec, tsim, mvec, hvec, nvec, isim2);

subplot(2,2,4);
plot(tsim, 1e3*vsim2, 'k');
hold on;
xlabel('Time (s)');
ylabel('Membrane Potental (mV)');
hold on;

% For getting train figure
% plot(tsim, vsim1, 'k');
% xlabel('Time (s)');
% ylabel('Membrane Potental (V)');

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
