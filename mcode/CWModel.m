classdef CWModel < handle
    %CWModel The dynamic model of a dynamic system governed by Clohessy-Wiltshire equation
    %   此处显示详细说明
    
    properties
        r = 5; % radius of the system
        p = [0; 0; 0]; % position of the system
        v = [0; 0; 0]; % velocity of the system, [p;v] renders the state of the system
        target = [0; 0; 0]; % target position
        history = [];
        stepNow = 0; % current step number
        isRecordHistory = true; % whether to record the history of the system
    end

    properties (Constant)
        % the angular velocity of the chasing satellite, say the chasing satellite is international space station
        % omega is then given by sqrt(mu/(R_earth+R_isshight)) = sqrt(3.986e14/(6.371+0.4)*1e6^3) = 0.00113
        % the period of ISS is then given by 2*pi/omega = 2*pi/(0.00113) = 5560s, which is nearly 93 minutes
        omega = 0.00113; 
    end

    % methods to manage class level properties (static properties), 
    % set dt with classname.setgetdt(<set value>), and get dt with classname.setgetdt()
    % see also: https://www.mathworks.com/help/releases/R2024a/matlab/matlab_oop/static-data.html#buvyy2x
    methods (Static)
        function out = setgetdt(dt_init)
           persistent dt;
           if nargin
              dt = dt_init;
           end
           out = dt;
        end

        function out = setgetsteps(steps_init)
           persistent steps;
           if nargin
              steps = steps_init;
           end
           out = steps;
        end
    end

    methods
        function obj = CWModel(p_init, v_init, target)
            %CWModel constructor
            obj.p = p_init;
            obj.v = v_init;
            obj.target = target;

            % if recored the history
            if obj.isRecordHistory
                obj.history = zeros(CWModel.setgetsteps(), 6);
                obj.recordOnce();
            end

            if isempty(CWModel.setgetdt())
                % set the default dt
                CWModel.setgetdt(0.1);
            end
            if isempty(CWModel.setgetsteps())
                % set the default steps
                CWModel.setgetsteps(2000);
            end
        end
        
        function obj = stepForward(obj, u)
            %stepForward update the state of the system, the dynamics of the model
            %   stepForward(obj, u) updates the state of the system with the given control input u
            %   we use Euler method to update the state of the system

            obj.p = obj.p + obj.v * CWModel.setgetdt();
            obj.v = obj.v + (obj.getCWf() + u) * CWModel.setgetdt();

            if obj.isRecordHistory
                % record the state of the system
                obj.recordOnce();
            end
        end

        function f = getCWf(obj)
            %getCWf get the drifting part of the dynamics of CW model
            %   the CW-model is given by
            %   x_dotdot = -2*omega*y_dot + u_x = f_x + u_x
            %   y_dotdot = 2*omega*x_dot + 3*omega^2*y + u_y = f_y + u_y
            %   z_dotdot = -omega^2*z + u_z = f_z + u_z
            % this function returns the drifting part of the dynamics
            %   f = [f_x; f_y; f_z]

            f = zeros(3, 1);
            f(1) = -2 * CWModel.omega * obj.v(2);
            f(2) = 2 * CWModel.omega * obj.v(1) + 3 * CWModel.omega^2 * obj.p(2);
            f(3) = -CWModel.omega^2 * obj.p(3);
        end

        function obj = recordOnce(obj)
            %recordOnce record the state of the system
            %   recordOnce(obj) records the state of the system
            %   in the history array.
            obj.stepNow = obj.stepNow + 1;
            obj.history(obj.stepNow, :) = [obj.p', obj.v'];
        end
    end

    methods (Static)
        function f = getAgentCWf(p, v)
            %getAgentCWf get the drifting part of the dynamics of CW model
            %   f = getAgentCWf(p, v) returns the drifting part of the dynamics
            %   of the system with the given position p and velocity v.
            %   f = [f_x; f_y; f_z]
            f = zeros(3, 1);
            f(1) = -2 * CWModel.omega * v(2);
            f(2) = 2 * CWModel.omega * v(1) + 3 * CWModel.omega^2 * p(2);
            f(3) = -CWModel.omega^2 * p(3);
        end
    end

    %% navigation control related

    methods 
        function PIDu = getPIDcontrol(obj)
            %getPIDcontrol get the control input of the system
            %   PIDu = getPIDcontrol(obj) returns the control input of the system
            %   based on the current state and target position.
            %   The control input is calculated using a PID controller.
            %   NOTE: the integral part is not implemented
            
            Kp = 0.004; % proportional gain
            Kd = 0.1; % derivative gain
            % Ki = 0.01; % integral gain
            
            e = obj.target - obj.p; % error
            de = -obj.v; % derivative of error
            
            PIDu = Kp * e + Kd * de;
        end

    end


    %% safe control related

    % ========================  distributed safety filtering =========================
    properties (SetAccess = private)
        % the safety control related parameters
        alpha1 = 0.05;
        alpha2 = 0.05;
        privArr = [];
    end

    methods
        function [u_safe, isQPfeasible] = distributedSafeFiltering(obj, u_ref, senseMat, privilege)
            %DISTRIBUTEDSAFEFILTERING get the safe-guarenteed control in a reciprocal way
            %   u_safe = reciprocalSafeFiltering(obj, u_ref, senseMat) returns the
            %   safe-guarenteed control input of the system based on the current state, 
            %   target position and the sensed data.
            %   u_ref is the reference control, which can be given by the PID controller.
            %   senseMat is a n*4 matrix, where n is the number of sensed objects, with each row
            %   the position and radius of the sensed object, [x, y, z, r].
            
            N_obj = size(senseMat, 1); % number of sensed objects
            obj.privArr = 0.5*ones(N_obj, 1); % the privilege array
            % set resArr to privilege if it is given
            if nargin > 3 && size(privilege, 1) == N_obj && size(privilege, 2) == 1
                obj.privArr = privilege;
            end
            
            [CBFconsA, CBFconsb] = obj.getCBFconstraints(senseMat); % get the CBF constraints
            % safety filtering, minimizing the distance to the reference control while satisfying the CBF constraints
            % u_safe = argmin ||u - u_ref||^2
            % s.t. CBFconsA * u <= CBFconsb
            H = 2*eye(3); % quadratic term
            f = -2*u_ref; % linear term
            [u_safe, modifValue, exitflag] = quadprog(H, f, CBFconsA, CBFconsb, [], [], [], [], [], CWModel.quadopt); % solve the quadratic programming problem

            % if the QP fails, use the minimum invasive result
            if exitflag < 0
                warning("QP error, using minimum invasive result")
                u_safe = CWModel.minmaxlinear(CBFconsA, CBFconsb); 
                isQPfeasible = 0;
            else
                isQPfeasible = 1;
            end
        end 


        function [CBFconsA, CBFconsb] = getCBFconstraints(obj, senseMat)
            %getCBFconstraints get the CBF constraints
            %   CBFcons = getCBFconstraints(obj, senseMat) returns the CBF constraints
            %   based on the current state and target position.
            %   senseMat is a n*7 matrix, where n is the number of sensed objects, with each row
            %   the position, velocity and radius of the sensed object, [x, y, z, vx, vy, vz, r].
            %   CBFconsA is a n*m matrix and CBFconsb is a n*1 matrix, where n is the number of sensed objects, m is the 
            %   demension of the control input.
            
            N_obj = size(senseMat, 1); % number of sensed objects
            CBFconsA = zeros(N_obj, 3); % CBF constraints
            CBFconsb = zeros(N_obj, 1); % CBF constraints
            
            for i = 1:N_obj
                p_i = senseMat(i, 1:3)';
                v_i = senseMat(i, 4:6)';
                f_i = CWModel.getAgentCWf(p_i, v_i); % get the drifting part of the dynamics of the sensed object
                r_i = senseMat(i, 7);
                R_i = r_i + obj.r;
                
                d_i = norm(obj.p - p_i); % distance between the system and the sensed object
                n_hat = (obj.p - p_i) / d_i; % unitvector pointing from the sensed object to the system
                v_rel = obj.v - v_i; % relative velocity
                
                % control of agent i can guarentee safety by satisfying CBF scaled by privilege
                CBFconsA(i, :) = -n_hat;
                
                % privilege matrix for only positive definite parts
                CBFconsb(i) = (obj.alpha1 + obj.alpha2) * n_hat' * obj.v  ...
                                + n_hat' * obj.getCWf() ...
                                + obj.privArr(i)* ( obj.alpha1*obj.alpha2* (d_i - R_i) ...
                                                    + 1/d_i * ( v_rel'*v_rel - n_hat'*v_rel * n_hat'*v_rel ) ...
                              );

                % privilege matrix for the whole b
                % CBFconsb(i) = obj.privArr(i)* ( (obj.alpha1 + obj.alpha2) * n_hat' * (obj.v - v_i)  ...
                %                                 +  obj.alpha1*obj.alpha2* (d_i - R_i) ...
                %                                 + 1/d_i * ( v_rel'*v_rel - n_hat'*v_rel * n_hat'*v_rel ) ...
                %                                 + n_hat' * (obj.getCWf() - f_i) ...
                %                                );

                if d_i < R_i - 1E-5
                    % give a warn of collision, collision may happen due to
                    % too large alpha and too large dt
                    warning("collision happened!");
                    d_i - R_i
                end
            end 
        end

        function [u_safe, isQPfeasible] = distributedSafeFilteringBenchmark(obj, u_ref, senseMat, privilege)
            %DISTRIBUTEDSAFEFILTERINGBENCHMARK get the safe-guarenteed control in a reciprocal way, assuming other agents are not controlled
            %   u_safe = reciprocalSafeFiltering(obj, u_ref, senseMat) returns the
            %   safe-guarenteed control input of the system based on the current state, 
            %   target position and the sensed data.
            %   u_ref is the reference control, which can be given by the PID controller.
            %   senseMat is a n*4 matrix, where n is the number of sensed objects, with each row
            %   the position and radius of the sensed object, [x, y, z, r].
            
            N_obj = size(senseMat, 1); % number of sensed objects
            obj.privArr = 0.5*ones(N_obj, 1); % the privilege array
            % set resArr to privilege if it is given
            if nargin > 3 && size(privilege, 1) == N_obj && size(privilege, 2) == 1
                obj.privArr = privilege;
            end
            
            [CBFconsA, CBFconsb] = obj.getCBFconstraintsBenchmark(senseMat); % get the CBF constraints
            % safety filtering, minimizing the distance to the reference control while satisfying the CBF constraints
            % u_safe = argmin ||u - u_ref||^2
            % s.t. CBFconsA * u <= CBFconsb
            H = 2*eye(3); % quadratic term
            f = -2*u_ref; % linear term
            [u_safe, modifValue, exitflag] = quadprog(H, f, CBFconsA, CBFconsb, [], [], [], [], [], CWModel.quadopt); % solve the quadratic programming problem

            % if the QP fails, use the minimum invasive result
            if exitflag < 0
                warning("QP error, using minimum invasive result")
                u_safe = CWModel.minmaxlinear(CBFconsA, CBFconsb); 
                isQPfeasible = 0;
            else
                isQPfeasible = 1;
            end
        end 

        function [CBFconsA, CBFconsb] = getCBFconstraintsBenchmark(obj, senseMat)
            %getCBFconstraints get the CBF constraints, assuming other agents are not controlled
            %   CBFcons = getCBFconstraints(obj, senseMat) returns the CBF constraints
            %   based on the current state and target position.
            %   senseMat is a n*7 matrix, where n is the number of sensed objects, with each row
            %   the position, velocity and radius of the sensed object, [x, y, z, vx, vy, vz, r].
            %   CBFconsA is a n*m matrix and CBFconsb is a n*1 matrix, where n is the number of sensed objects, m is the 
            %   demension of the control input.
            
            N_obj = size(senseMat, 1); % number of sensed objects
            CBFconsA = zeros(N_obj, 3); % CBF constraints
            CBFconsb = zeros(N_obj, 1); % CBF constraints
            
            for i = 1:N_obj
                p_i = senseMat(i, 1:3)';
                v_i = senseMat(i, 4:6)';
                f_i = CWModel.getAgentCWf(p_i, v_i); % get the drifting part of the dynamics of the sensed object
                r_i = senseMat(i, 7);
                R_i = r_i + obj.r;
                
                d_i = norm(obj.p - p_i); % distance between the system and the sensed object
                n_hat = (obj.p - p_i) / d_i; % unitvector pointing from the sensed object to the system
                v_rel = obj.v - v_i; % relative velocity
                
                % control of agent i can guarentee safety by satisfying CBF scaled by privilege
                CBFconsA(i, :) = -n_hat;
            
                CBFconsb(i) = (obj.alpha1 + obj.alpha2) * n_hat' * (obj.v - v_i)  ...
                                + n_hat' * (obj.getCWf() - f_i) ...
                                + obj.alpha1*obj.alpha2* (d_i - R_i) ...
                                + 1/d_i * ( v_rel'*v_rel - n_hat'*v_rel * n_hat'*v_rel );

                if d_i < R_i - 1E-5
                    % give a warn of collision, collision may happen due to
                    % too large alpha and too large dt
                    warning("collision happened!");
                    d_i - R_i
                end
            end 
        end

    end

    % ======================== centralized safety filtering ========================
    properties (Constant)
        % the safety control related parameters
        % class kappa function for centralized safety filtering,
        % should be the same as alpha1 and alpha2 of class instances
        alpha1_c = 0.05;
        alpha2_c = 0.05;
    end

    methods (Static)
        function U_safe = centralizedSafeFiltering(U_ref, globalSenseMat)
            %CENTRALIZEDSAFEFILTERING get the safe-guaranteed control in a centralized way
            %   U_safe = centralizedSafeFiltering(U_ref, globalSenseMat) does safety filtering 
            %   simultaneously for all agents
            %   U_ref is the concatenated reference control as U_ref = [u_ref1; u_ref2; ...; u_refN]
            %   globalSenseMat contains all position, velocity and radius data of all agents, managed like 
            %   [x1, y1, z1, vx1, vy1, vz1, r1; x2, y2, z2, vx2, vy2, vz2, r2; ...; xN, yN, zN, vxN, vyN, vzN, rN] in R^(N*7)
            %   U_safe is the concatenated safe-guaranteed control as U_safe = [u_safe1; u_safe2; ...; u_safeN] in R^(3*N)
            %   No privillage needed since the optimization is done in a global way

            N = size(globalSenseMat, 1); % number of agents
            CBFconsA = zeros(0.5*N*(N-1), 3*N); % CBF constraints, N(N-1)/2 constriants should be considered
            CBFconsb = zeros(0.5*N*(N-1), 1); % CBF constraints, right hand side

            % get the CBF constraints
            for i = 1:N
                % get the sensed data of agent i
                p_i = globalSenseMat(i, 1:3)';
                v_i = globalSenseMat(i, 4:6)';
                f_i = CWModel.getAgentCWf(p_i, v_i); % get the drifting part of the dynamics of the sensed object
                r_i = globalSenseMat(i, 7);

                for j = i+1:N
                    % get the sensed data of agent j
                    p_j = globalSenseMat(j, 1:3)';
                    v_j = globalSenseMat(j, 4:6)';
                    f_j = CWModel.getAgentCWf(p_j, v_j); % get the drifting part of the dynamics of the sensed object
                    r_j = globalSenseMat(j, 7);
                    R = r_i + r_j;
                    
                    % get the distance between the two agents
                    d_ij = norm(p_i - p_j); % distance between the two agents
                    n_hat_ij = (p_i - p_j) / d_ij; % unitvector pointing from agent j to agent i
                    v_ij = v_i - v_j; % relative velocity
                    
                    % CBF constraints
                    % the row index is given by the sum of 
                    % 1. number of constraints of 1 to i-1 agent, that is (N-1) + (N-2) + ... + (N-(i-1)) = (i-1)*(2*N-i)/2 = N(i-1) - (i-1)*(i)/2
                    % 2. number of constraints from i+1 to j, that is j-i
                    CBFconsA(0.5*(i-1)*(2*N-i) + j-i, (i-1)*3+1:(i-1)*3+3) = -n_hat_ij; % for agent i and block i
                    CBFconsA(0.5*(i-1)*(2*N-i) + j-i, (j-1)*3+1:(j-1)*3+3) = n_hat_ij; % for agent j and block j
                    % right hand side
                    CBFconsb(0.5*(i-1)*(2*N-i) + j-i) = (CWModel.alpha1_c + CWModel.alpha2_c) * n_hat_ij' * v_ij ...
                                                        + CWModel.alpha1_c*CWModel.alpha2_c* (d_ij - R) ...
                                                        + 1/d_ij * (v_ij'*v_ij - n_hat_ij'*v_ij * n_hat_ij'*v_ij) ...
                                                        + n_hat_ij' * (f_i - f_j);

                    if d_ij < R - 1E-5
                        % give a warn of collision, collision may happen due to
                        % too large alpha and too large dt
                        warning("collision happened!");
                        d_ij - R
                    end
                end
            end

            % do CBF-QP as
            % U_safe = argmin ||U - U_ref||^2
            % s.t. CBFconsA * U <= CBFconsb
            H = 2*eye(3*N); % quadratic term
            f = -2*U_ref; % linear term
            [U_safe, modifValue, exitflag] = quadprog(H, f, CBFconsA, CBFconsb, [], [], [], [], [], CWModel.quadopt); % solve the quadratic programming problem

            % if the QP fails, use the minimum invasive result
            if exitflag < 0
                warning("QP error, using minimum invasive result")
                U_safe = CWModel.minmaxlinear(CBFconsA, CBFconsb); 
            end
        end
    end

    %% under testing part: privilege optimizing
    methods (Static)
        function newPrivilegeMat = privilegeOptimze1(U_safe, SenseMat, privilegeMat)
            %PRIVILEGEOPTMIZE optimize the privilege matrix with the centralized safety control
            %   newPrivilegeMat = privilegeOptimze(U_safe, SenseMat, privilegeMat) returns the optimized
            %   privilege matrix by making the decentrialized safety constraints satisfied constrained by the
            %   centralized safety output.
            %   U_safe is the concatenated safe-guaranteed control as U_safe = [u_safe1; u_safe2; ...; u_safeN] in R^(3*N)
            %   SenseMat is a n*7 matrix, where n is the number of sensed objects, with each row
            %   the position, velocity and radius of the sensed object, [x, y, z, vx, vy, vz, r].
            %   privilegeMat is the privilege matrix, which is a n*n matrix, where n is the number of agents,
            %   we define (p_ij + p_ji) = 1, p_ij in [0, 1]. The larger p_ij, the more easy the local safety constraint to be satisfied.
    
            N = size(SenseMat, 1); % number of agents
            newPrivilegeMat = privilegeMat; % initialize the new privilege matrix
            
            for i = 1:N
                for j = i+1:N
                    % get the sensed data of agent i
                    p_i = SenseMat(i, 1:3)';
                    v_i = SenseMat(i, 4:6)';
                    r_i = SenseMat(i, 7);
                    u_sci = U_safe((i-1)*3+1:(i-1)*3+3);
    
                    % get the sensed data of agent j
                    p_j = SenseMat(j, 1:3)';
                    v_j = SenseMat(j, 4:6)';
                    r_j = SenseMat(j, 7);
                    u_scj = U_safe((j-1)*3+1:(j-1)*3+3);
                    
                    R = r_i + r_j;
                    
                    d_ij = norm(p_i - p_j); % distance between the two agents
                    n_hat_ij = (p_i - p_j) / d_ij; % unitvector pointing from agent j to agent i
                    n_hat_ji = -n_hat_ij;
                    v_ij = v_i - v_j; % relative velocity
                    v_ji = -v_ij;
                    
                    % TODO: do QP to optimize the privilege matrix
                    Acons = zeros(4, 1);
                    bcons = zeros(4, 1);
                    Acons(1) =  - ( CWModel.alpha1_c*CWModel.alpha2_c * (d_ij - R) + 1/d_ij * (v_ij'*v_ij - n_hat_ij'*v_ij * n_hat_ij'*v_ij) );
                    Acons(2) = ( CWModel.alpha1_c*CWModel.alpha2_c * (d_ij - R) + 1/d_ij * (v_ij'*v_ij - n_hat_ij'*v_ij * n_hat_ij'*v_ij) );
                    Acons(3) = -1;
                    Acons(4) = 1;
                    bcons(1) = n_hat_ij'*u_sci + (CWModel.alpha1_c + CWModel.alpha2_c) * n_hat_ij'*v_i;
                    bcons(2) = n_hat_ji'*u_scj + (CWModel.alpha1_c + CWModel.alpha2_c) * n_hat_ji'*v_j + CWModel.alpha1_c*CWModel.alpha2_c * (d_ij - R) + 1/d_ij * (v_ij'*v_ij - n_hat_ij'*v_ij * n_hat_ij'*v_ij);
                    bcons(3) = 0;
                    bcons(4) = 1;
    
                    % solve the QP problem
                    H = 2*eye(1); % quadratic term
                    f = -2*privilegeMat(i, j); % linear term
                    [p_ij, modifValue, exitflag] = quadprog(H, f, Acons, bcons, [], [], [], [], [], CWModel.quadopt); % solve the quadratic programming problem
                    
                    if exitflag >= 0
                        newPrivilegeMat(i,j) = p_ij;
                        newPrivilegeMat(j,i) = 1-p_ij;
                    end
                end
            end
        end

        function newPrivilegeMat = privilegeOptimize2(senseMat, privilegeMat, Goals)
            %PRIVILEGEOPTMIZE optimize the privilege matrix to minimize the difference of safety output between the centralized and decentralized safety filter
            %   newPrivilegeMat = privilegeOptimize(senseMat, privilegeMat) returns the optimized
            %   privilege matrix by tuning the privilege matrix to minimize the difference of safety output between the centralized and decentralized safety filter
            %   that is min_P ||U_safe_decentralized - U_safe_centralized||^2, with U_safe_decentralized and U_safe_centralized given by
            %   decentralizedSafeFiltering and centralizedSafeFiltering method, respectively.
            %   senseMat is a n*7 matrix, where n is the number of sensed objects, with each row
            %   the position, velocity and radius of the sensed object, [x, y, z, vx, vy, vz, r].   
            %   privilegeMat is the privilege matrix, which is a n*n matrix, where n is the number of agents,
            %   we define p_ij in [0, 1] and p_ij + p_ji <= 1, p_ii = 0. The larger p_ij, the more easy the local safety constraint to be satisfied.

            N = size(senseMat, 1); % number of agents
            agentCell = cell(N, 1); % initialize the agent cell
            U_ref = zeros(N*3, 1); % initialize the reference control
            for i = 1:N
                % get the sensed data of agent i
                p_i = senseMat(i, 1:3)';
                v_i = senseMat(i, 4:6)';
                agentCell{i} = CWModel(p_i, v_i, Goals(:, i)); % create the agent cell
                u_ref = agentCell{i}.getPIDcontrol();
                U_ref((i-1)*3+1:(i-1)*3+3) = u_ref; % get the reference control
            end
            U_safe_c = CWModel.centralizedSafeFiltering(U_ref, senseMat);

            privilegeMatArr = privilegeMat2Arr(privilegeMat); % convert the privilege matrix to array, so the optimization has no constraints

            % optimized the privilege matrix
            % newPrivilegeMat = fmincon(@getU_diff, privilegeMat, [], [], [], [], zeros(N, N), ones(N, N), @nonlcon, CWModel.fminconopt);
            newPrivilegeMatArr = fminunc(@getU_diff, privilegeMatArr, CWModel.fminuncopt); % optimize the privilege matrix
            % newPrivilegeMatArr = fminsearch(@getU_diff, privilegeMatArr, CWModel.fminsearchopt); % optimize the privilege matrix with gradient-free method

            newPrivilegeMat = arr2privilegeMat(newPrivilegeMatArr); % convert the array to privilege matrix

            function U_diff = getU_diff(PrivMatArr)
                % get the difference of safety output with the given privilege matrix as variable
                % convert the array to privilege matrix
                PrivMat = arr2privilegeMat(PrivMatArr);
                [U_safe_d, n_notfeas] = getDecentralizedU(PrivMat);
                U_diff = norm(U_safe_d - U_safe_c) + 1E5*n_notfeas; % add penalty for not feasible case
            end

            function [U_safe_d, n_notfeas] = getDecentralizedU(PrivMat)
                % get the decentralized safety output with the given privilege matrix as variable
                U_safe_d = zeros(N*3, 1);
                n_notfeas = 0;
                for j = 1:N
                    senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                    u = U_ref((j-1)*3+1:(j-1)*3+3); % reference control
                    privilegeArray = PrivMat(j, :);
                    privilegeArray(j) = [];
                    [u_safe, isQPfeasible] = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray');
                    U_safe_d((j-1)*3+1:(j-1)*3+3) = u_safe;
                    if ~isQPfeasible
                        n_notfeas = n_notfeas + 1;
                    end
                end
            end

            function [c, ceq] = nonlcon(PrivMat)
                % non-linear constraints for the privilege matrix
                % c = PrivMat' + PrivMat - 1; % p_ij + p_ji <= 1
                % ceq = diag(PrivMat); % p_ii = 0

                c = [];
                ctemp = PrivMat' + PrivMat - (ones(N) - eye(N)); % p_ij + p_ji = 1
                ceq = [diag(PrivMat);
                       ctemp(:)]; % p_ii = 0
            end

            function privilegeMatArr = privilegeMat2Arr(PrivMat)
                % convert the privilege matrix to array
                privilegeMatArr = zeros(N*(N-1)/2, 1);
                index = 1;
                for ii = 1:N
                    for jj = ii+1:N
                        privilegeMatArr(index) = PrivMat(ii, jj);
                        index = index + 1;
                    end
                end
            end

            function PrivMat = arr2privilegeMat(privilegeMatArr)
                % convert the array to privilege matrix
                PrivMat = zeros(N, N);
                index = 1;
                for ii = 1:N
                    for jj = ii+1:N
                        PrivMat(ii, jj) = privilegeMatArr(index);
                        PrivMat(jj, ii) = 1 - privilegeMatArr(index);
                        index = index + 1;
                    end
                end
            end
        end


        function newPrivilegeMat = privilegeOptimize(senseMat, privilegeMat, Goals)
            %PRIVILEGEOPTMIZE optimize the privilege matrix to minimize the difference between reference control U_ref and distributed control U_safe_d, that is
            %   min_P ||U_safe_decentralized - U_ref||^2
            %   senseMat is a n*7 matrix, where n is the number of sensed objects, with each row
            %   the position, velocity and radius of the sensed object, [x, y, z, vx, vy, vz, r].   
            %   privilegeMat is the privilege matrix, which is a n*n matrix, where n is the number of agents,
            %   we define p_ij in [0, 1] and p_ij + p_ji == 1, p_ii = 0. The larger p_ij, the more easy the local safety constraint to be satisfied.

            N = size(senseMat, 1); % number of agents
            agentCell = cell(N, 1); % initialize the agent cell
            U_ref = zeros(N*3, 1); % initialize the reference control
            for i = 1:N
                % get the sensed data of agent i
                p_i = senseMat(i, 1:3)';
                v_i = senseMat(i, 4:6)';
                agentCell{i} = CWModel(p_i, v_i, Goals(:, i)); % create the agent cell
                u_ref = agentCell{i}.getPIDcontrol();
                U_ref((i-1)*3+1:(i-1)*3+3) = u_ref; % get the reference control
            end

            privilegeMatArr = privilegeMat2Arr(privilegeMat); % convert the privilege matrix to array, so the optimization has no constraints

            % optimized the privilege matrix
            % newPrivilegeMat = fmincon(@getU_diff, privilegeMat, [], [], [], [], zeros(N, N), ones(N, N), @nonlcon, CWModel.fminconopt);
            newPrivilegeMatArr = fminunc(@getU_diff, privilegeMatArr, CWModel.fminuncopt); % optimize the privilege matrix

            newPrivilegeMat = arr2privilegeMat(newPrivilegeMatArr); % convert the array to privilege matrix

            function U_diff = getU_diff(PrivMatArr)
                % convert the array to privilege matrix
                PrivMat = arr2privilegeMat(PrivMatArr);
                % get the difference of safety output with the given privilege matrix as variable
                U_safe_d = getDecentralizedU(PrivMat);
                U_diff = norm(U_safe_d - U_ref);
            end

            function U_safe_d = getDecentralizedU(PrivMat)
                % get the decentralized safety output with the given privilege matrix as variable
                U_safe_d = zeros(N*3, 1);
                for j = 1:N
                    senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                    u = U_ref((j-1)*3+1:(j-1)*3+3); % reference control
                    privilegeArray = PrivMat(j, :);
                    privilegeArray(j) = [];
                    u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray');
                    U_safe_d((j-1)*3+1:(j-1)*3+3) = u_safe;
                end
            end

            function [c, ceq] = nonlcon(PrivMat)
                % non-linear constraints for the privilege matrix
                % c = PrivMat' + PrivMat - 1; % p_ij + p_ji <= 1
                % ceq = diag(PrivMat); % p_ii = 0
                
                c = [];
                ctemp = PrivMat' + PrivMat - (ones(N) - eye(N)); % p_ij + p_ji <= 1
                ceq = [diag(PrivMat);
                       ctemp(:)]; % p_ii = 0
            end

            function privilegeMatArr = privilegeMat2Arr(PrivMat)
                % convert the privilege matrix to array
                privilegeMatArr = zeros(N*(N-1)/2, 1);
                index = 1;
                for ii = 1:N
                    for jj = ii+1:N
                        privilegeMatArr(index) = PrivMat(ii, jj);
                        index = index + 1;
                    end
                end
            end

            function PrivMat = arr2privilegeMat(privilegeMatArr)
                % convert the array to privilege matrix
                PrivMat = zeros(N, N);
                index = 1;
                for ii = 1:N
                    for jj = ii+1:N
                        PrivMat(ii, jj) = privilegeMatArr(index);
                        PrivMat(jj, ii) = 1 - privilegeMatArr(index);
                        index = index + 1;
                    end
                end
            end
        end

        
    end


    %% optimization related
    properties (Constant)
        quadopt = optimoptions("quadprog", "Display", "off");
        fminconopt = optimoptions("fmincon", "Display", "off");
        fminuncopt = optimoptions("fminunc", "Display", "off");
        fminsearchopt = optimset('Display','off');
    end

    methods (Static)
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
    end


    %% visualization related

    methods 
        function plotHistory(obj, varargin)
            %plotHistory plot the history of the system with 3D quiver arrows
            %   plotHistory(obj, stepInterval) plots the trajectory of the system
            %   in 3D space with velocity vectors. stepInterval specifies the
            %   interval at which arrows are plotted (default is 30).
    
            % Define default values
            defaultStepInterval = 3;
            defaultColor = [];

            % Set up input parser
            ip = inputParser;
            addParameter(ip, 'stepInterval', defaultStepInterval, @(x) isscalar(x) && x > 0);
            addParameter(ip, 'color', defaultColor, @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
            parse(ip, varargin{:});

            % Extract parameters
            stepInterval = ip.Results.stepInterval;
            color = ip.Results.color;
    
            % Plot the trajectory
            p = plot3(obj.history(:, 1), obj.history(:, 2), obj.history(:, 3), 'LineWidth', 1.5, 'LineStyle', 'none');

            if ~isempty(color)
                p.Color = color; % Set the color of the trajectory
            end
    
            % Prepare data for quiver arrows
            xs = obj.history(1:stepInterval:end, 1);
            ys = obj.history(1:stepInterval:end, 2);
            zs = obj.history(1:stepInterval:end, 3);
            vxs = obj.history(1:stepInterval:end, 4);
            vys = obj.history(1:stepInterval:end, 5);
            vzs = obj.history(1:stepInterval:end, 6);

            % plot the trajectory with scatters
            s = scatter3(xs, ys, zs);
            s.Marker = "o";
            s.SizeData = 15;
            s.MarkerFaceColor = "flat";
            s.MarkerEdgeColor = "flat";
            s.AlphaDataMapping = "scaled";
            s.MarkerFaceColor = p.Color;
            s.MarkerEdgeColor = p.Color;
            s.MarkerFaceAlpha = 0.20;
            s.MarkerEdgeAlpha = 0.20;
   
            quiverMultiplier = 12;
            vnorm = (vxs.^2 + vys.^2 + vzs.^2).^(0.5);
            vxs = vxs./vnorm * quiverMultiplier;
            vys = vys./vnorm * quiverMultiplier;
            vzs = vzs./vnorm * quiverMultiplier;
    
            % Plot velocity vectors using quiver3
            quiver3(xs(1:2), ys(1:2), zs(1:2), vxs(1:2), vys(1:2), vzs(1:2), 'filled', 'AutoScale', 'off', 'Color', p.Color, 'LineWidth', 2.0, 'MaxHeadSize', 30, 'LineStyle','-');
    
            % Highlight the target position if goal is arrived
            plot3(obj.target(1), obj.target(2), obj.target(3), 'Marker', 'pentagram',  'MarkerSize', 10, 'LineWidth', 2, 'Color', p.Color, 'LineStyle', 'none');
            % plot3(obj.history(1,1), obj.history(1,2), obj.history(1,3)-5, 'Marker', '^',  'MarkerSize', 5, 'LineWidth', 2, 'Color', p.Color);
   
        end
    end
end

