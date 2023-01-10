
function [ i1vec, i2vec, tvec ] = generate_iapp(started, ended, dt, mu, S_bar, tau_c)
% returns 2 randomly generated Iapp vector, each fed into neuron 1 and 2.
% the time windows (period vector) for constant Iapp should be chosen 
% from an exponential distribution, with a mean of tau_c (30ms default)
% sum of all periods should be equal to tmax. 
% the returned ivec is "voltage term, since resistance is multiplied.
% when incorporating in main code for H-H model, divide by resistance term (Rm)

% initiate tvec
tvec = started:dt:ended;                                
tmax = ended-started;                   % tmax indicates total running time

% Make period vector (time windows)
period_vec = zeros(1, length(tvec));    % lengthy initialization (size t)
flag = false;
i = 1;
sum = 0;

while ~flag
    period = exprnd(tau_c);
    sum = sum + period;
    if sum <= tmax
        period_vec(i) = period;
        i = i + 1;
    else
        flag = true;
    end
end


period_vecy = transpose(nonzeros(period_vec));           % remove all 0. only valid values. transpose to row vector

% Fill last element of period_vec
track = 0;
for k = 1:length(period_vecy)
    track = track + period_vecy(k);
end

period_vecy(end+1) = tmax-track;                        % fill last period to match tmax
cum_period = cumsum(period_vecy);                       % cumulative vector of periods. last element = tmax

% Make p and s vector
p1vec = zeros(1, length(tvec));
p2vec = zeros(1, length(tvec));
s1vec = zeros(1, length(tvec));
s2vec = zeros(1, length(tvec));

% Fill p and s vector
for k = 1:length(cum_period)-1
    rand_p1 = rand*S_bar;
    rand_p2 = rand*S_bar;
    rand_s1 = rand*S_bar;
    
    for l = (floor(cum_period(k)/dt)+1):(floor(cum_period(k+1)/dt))
        p1vec(l) = rand_p1;
        p2vec(l) = rand_p2;
        s1vec(l) = rand_s1;
        s2vec(l) = S_bar - s1vec(l);
    end
end

rand_py1 = rand*S_bar;
rand_py2 = rand*S_bar;
rand_sy1 = rand*S_bar;

for l = 1:floor(cum_period(1)/dt)           % fill start~first period
    p1vec(l) = rand_py1;
    p2vec(l) = rand_py2;
    s1vec(l) = rand_sy1;
    s2vec(l) = S_bar - s1vec(l);
end

% Make I vector, linear combination of P and S
i1vec = zeros(1, length(tvec));
i2vec = zeros(1, length(tvec));

for i = 1:length(tvec)
    i1vec(i) = ((1-mu)*p1vec(i)) + (mu*s1vec(i));
    i2vec(i) = ((1-mu)*p2vec(i)) + (mu*s2vec(i));
end
end

