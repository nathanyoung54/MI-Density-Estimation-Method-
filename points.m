
function [ result ] = points(matrix, i, n, h)
% i is standard point
% n is length of matrix
% for a standard point i, generate 'h' close points and return their
% indices

lhs = {};
rhs = {};

for j=1:i-1
    save1 = [matrix(i, j)];
    lhs{end+1} = save1;
end

for j =i+1:n
    save2 = [matrix(i, j)];
    rhs{end+1} = save2;
end

excised = cell2mat([lhs rhs]);
excised_sorted = sort(excised);
excised_sorted = excised_sorted(1:h-1);         % find h-1 close points, 1 is itself

result = {};

while length(result) < h-1
    for value=1:length(excised_sorted)
    index = find(excised==excised_sorted(value));
    result{end+1} = index;
    end
end

result = cell2mat(result);
result = result(1:h-1);

end