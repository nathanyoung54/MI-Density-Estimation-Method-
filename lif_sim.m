
function [v_sim, find_spike_sim] = lif_sim(t, ref_period, I_i)
% LIF (Leaky-integrate-and-fire) model
% % returns simulated membrane potential vector and spike train vector
% spike train vector indicates "time points" where spikes occured
% t is time vector
% ref_period is fixed time for no spikes
% I_i is voltage term, Rm*Iapp, generated from generate_iapp

global EL Vth Vre tau_m
global dt_lif 

v_sim = zeros(size(t));                                     % MP vector same size with time vector
v_sim(1) = EL;                                              % initialize with leak potential
spike_sim = zeros(size(t));                                 % save spike-occuring t

time_record = zeros(size(t));                               % create vector for time recording
time_record(1) = -ref_period*1.01;                          % initialize less than -ref_period

% Run simulation

for i = 2:length(t)                                         % loop through time points                       
    tval = t(i);                                            % value of time point
    if tval < time_record(i-1) + ref_period                 % tval needs to exceed refractory period
        v_sim(i) = Vre;
        time_record(i) = time_record(i-1);
    else
        dvdt = (EL - v_sim(i-1) + I_i(i))/tau_m;   
        v_sim(i) = v_sim(i-1) + dvdt*dt_lif;
        if v_sim(i) > Vth                                   % if MP > threshold, set to Vre
            spike_sim(i) = i*dt_lif;                        % save spike-occuring t
            v_sim(i) = Vre;
            time_record(i) = tval;                          % record the time for Vm > Vth
            v_sim(i-1) = 50e-3;                             % spike membrane potential 
        end                                                 % if tval does not exceed refractory period, Vm = Vre
    end
end

find_spike_sim = dt_lif*find(spike_sim);

end