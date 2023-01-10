
function [ info ] = information_from_matrix(u_matrix, v_matrix, u_h, v_h, leave_out)
% returns calculated mutual information from two spike train distance matrices
% u_h and v_h are smoothing factors (u_h = v_h)
% leave_out is set 1 (adding itself when considering neighbors)

import distance_matrix.* points.*

if ~exist('leave_out', 'var'), leave_out = 1; end

n = length(u_matrix);            % nxn matrix (u and v should have same length)
h = [u_h v_h];                   % to store h

u_size = size(u_matrix);
v_size = size(v_matrix);

% if size is not correct
if n ~= u_size(2) | n ~= v_size(1) | n ~= v_size(2)
    disp("wrong input!");
    exit;
end

% initiate variables
info = 0;
total_added = 0;

% calculate MI                  
for i=1:n
    u_points = points(u_matrix, i, n, u_h);
    v_points = points(v_matrix, i, n, v_h);
    
    % find intersection between u_points and v_points
    inter = intersect(u_points, v_points);
    hash_inter = leave_out + length(inter);
    
    if hash_inter>0
        info = info + log2((n+leave_out-1)*hash_inter/(u_h*v_h));
        total_added = total_added + 1;
    end

end

info = info/total_added;

end