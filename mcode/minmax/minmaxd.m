function optx = minmaxd(A, b)
% optimize the min max di(x) problem, where di is the distance (ESDF
% actually) to the planes A*x <= b
% this function is tested by testminmaxd

opt = optimoptions("linprog", "Display","off", "ConstraintTolerance", 1E-10, "OptimalityTolerance", 1E-10, "MaxIterations", 1E6, "MaxTime", 1E3);
EPS = eps;
M = 1E5;

% the norm of each row of A
aNorms = sqrt(sum(A.^2, 2));

% normalize, necessary since the direction vector should be of same norm
A_normd = A ./ aNorms;
b_normd = b ./ aNorms;

% add a bounding box to avoid unbounded
% -M <= x_i <= M
dim = size(A,2);
A_normd = [eye(dim);
          -eye(dim);
           A_normd];
b_normd = [ones(dim,1) * M;
           ones(dim,1) * M;
           b_normd];

% initial d
d = -M;
% initial point need to be on the line
x_star = ones(dim,1) * 0;
% x_star = sqrt(dim) * M * ones(dim,1);
% x_star = rand(dim,1);

% checkFeasibleFlag = true;
for i = (2*dim + 1):size(A_normd,1)
    % TODO: NEED CHECK THIS, FIND WHY INOPTIMAL EXISTS
    % always continue if constraints is feasible
    % if checkFeasibleFlag
    %     % solve LP
    %     [lp_x_star, ~, flag] = linprog(ones(dim,1), A_normd(1:i,:), b_normd(1:i,:), [], [], [], [], opt);
    %     if flag ~= -2 % not infeasible
    %         if flag == 1 % if converged
    %             x_star = lp_x_star; % let x_star on the boundary
    %         end
    %         continue
    %     else % == -2, infeasible
    %         checkFeasibleFlag = false; % not need to check feasibility anymore
    %     end
    % end
    if A_normd(i,:) * x_star - b_normd(i) < d
        continue
    end
    projA = [];
    projb = [];
    for j = 1:i-1
        % generate the equal distance line
        if abs(det([A_normd(i,:); A_normd(j,:)])) < EPS % almost parallel
            if A_normd(i,:) * A_normd(j,:)' > 0
                continue;
            else
                % checked, should use "-" to get c
                a = (A_normd(j,:) - A_normd(i,:)) / norm(A_normd(j,:) - A_normd(i,:)); 
                c = (b_normd(j) - b_normd(i)) / 2;
            end
        else
            a = (A_normd(j,:) - A_normd(i,:)) / norm(A_normd(j,:) - A_normd(i,:));
            c = (b_normd(j) - b_normd(i)) / norm(A_normd(j,:) - A_normd(i,:));
        end
        projA = [projA; a];
        projb = [projb; c];
    end
    % PROBLEM: THIS MAY BE UNBOUNDED, use bounding box here
    [x_star, ~, flag2] = linprog(A_normd(i,:), projA, projb, [], [], [], [], opt);
    d = A_normd(i,:) * x_star - b_normd(i);
    if flag2 ~= 1
        warning("LP not converged, may lead to error")
    end
end
optx = x_star;
end