
function [ fragments ] = chop_train(train, train_length, window_length, start_time, en_time);
% chops the spike train into size of "window_length"
% returns all the fragments in a cell array "fragments"

next = 1;
time = start_time;
fragments = {};

while time < en_time
    this_fragment = [];
    while (next <= train_length) & (train(next) < time + window_length)
        this_fragment(end+1) = train(next)-time;
        next = next + 1;
    end
    fragments{end+1} = this_fragment;
    time = time + window_length;
end
end