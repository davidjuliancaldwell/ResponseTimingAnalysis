X = rand(2,3,4,5,6);
nv   = ndims(X);
v    = ones(1, nv);
vLim = size(X);
ready = false;
while ~ready
    % Do what you need with X and the index vector v
    ...
        % Update the index vector:
    ready = true;       % Assume that the WHILE loop is ready
    for k = 1:nv
        v(k) = v(k) + 1;
        if v(k) <= vLim(k)
            ready = false;  % No, WHILE loop is not ready now
            break;          % v(k) increased successfully, leave "for k" loop
        end
        v(k) = 1;         % Reset v(k), proceed to next k
    end
end

3