function [ h, old_h ] = golden_h(mat1, mat2, old_h, biggest_h)
% Golden Mean Search for h

import information_from_matrix.* background.* 

n1 = length(mat1);
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
old_h = h;