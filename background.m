
function [ info ] = background(n, h)
% returns the calculated bias that needs to be subtracted from MI

import calculate_coefficient.*

info = 0;

for r=1:h
    prob = nchoosek(h-1,r-1) * nchoosek(n-h,h-r) / nchoosek(n-1,h-1);
    info = info + prob*log2(n*r/h^2);

end
end

% info = info + calculate_coefficient(n-1,h-1,r)*log2(n*r/(h*h));
%     prob = (factorial(h-1)/factorial(r-1)*factorial(h-r)) ...
%         * (factorial(n-h)/factorial(h-r)*factorial(n+r))...
%         / (factorial(n-1)/factorial(h-1)*factorial(n-h));
%     info = info + prob*log2(n*r/h*h);