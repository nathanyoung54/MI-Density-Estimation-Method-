
function [ tvec, ivec ] = pulsemaker(istart, ibase, npulses,...
    tau_sep, tau_int, tau_len, sigma_sep, sigma_int, sigma_len, tmax, dt)
% Pulse input current generator for Experiment #2,3 and 4
% istart: Default time applied current starts
% ibase: Default baseline current before/after pulse
% npulses: Default number of current pulses
% intensity: Default intensity
% tau_sep: mean of interval between pulses
% tau_int: mean of pulse intensity 
% tau_len: mean of one pulse length
% sigma_sep: sd of interval between pulses
% sigma_sep: sd of interval between pulses
% sigma_sep: sd of interval between pulses

% Initiate vectors
tvec = 0:dt:tmax;
ivec = ibase*ones(size(tvec));       % Initialize current vector at baseline

for pulse=1:npulses
    if sigma_sep==0
        pulsesep = tau_sep;
    else
        pulsesep = normrnd(tau_sep, sigma_sep);
    end
    if sigma_int==0
        intensity = tau_int;
    else
        intensity = normrnd(tau_int, sigma_int);
    end
    if sigma_len==0
        ilength = tau_len;
    else
        ilength = normrnd(tau_len, sigma_len);
    end
    pulsestart = istart + (pulse-1)*pulsesep;   % Onset time of pulse
    pulsestop = pulsestart + ilength;           % Offset time of pulse
    
    % make applied current a new value for duration of current pulse
    for q=round(pulsestart/dt):round(pulsestop/dt)
        ivec(q) = intensity;
    end
end