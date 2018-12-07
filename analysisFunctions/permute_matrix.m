function output = permute_matrix(input,dimInt)
% function to permute matrix along one of the dimensions of a 3D data
% matrix
% input -
% David.J.Caldwell
lengthShuffle = size(input,dimInt);

operateOn = [2,3];  % Operate on submatrix along 3rd and 5th dimension
v   = 1:ndims(X);
XX  = permute(X, [operateOn, setdiff(v, operateOn)]);
sXX = size(XX);
XX  = reshape(XX, [sXX(1), prod(sXX(2:end))]);
for k = 1:size(XX, 3)
    submatrix = XX(:, :, k);
    ... do what you want
end

end