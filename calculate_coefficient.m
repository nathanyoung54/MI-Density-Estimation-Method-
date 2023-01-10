function [ value ] = calculate_coefficient(n, h, r)
value = 1.0;

if r == 0
    for counter=0:h-1
        value = value * (n-h-counter)/(n-counter);
    end
end

for counter=0:h-r-1
    value = value * ((n-h-counter)/(n-counter))*((h-counter)/(h-r-counter));
end

for counter=0:r-1
    value = value * ((h-counter)/(n-h+r-counter));
end

