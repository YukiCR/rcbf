function x_star = minmaxlinear(A, b)
% optimize the min max di(x) problem, where di is the distance (ESDF
% actually) to the planes A*x <= b
% this function is tested by testminmaxd and faster than MATLAB fminimax

% the norm of each row of A
aNorms = sqrt(sum(A.^2, 2));

[n,dim] = size(A);
A_extend = [A, -ones(n, 1) .* aNorms];
b_extend = b;

opt = optimoptions("linprog", "Display","off");

f = zeros(dim+1, 1);
f(end) = 1;
[optx, ~, flag] = linprog(f, A_extend, b_extend, [], [], [], [], opt);
if flag ~= 1
    warning("LP not converged, may lead to error")
end

x_star = optx(1:end-1);

end