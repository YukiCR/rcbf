classdef twoOrderModel < handle
    %TWOORDERMODEL The dynamic model of a two-order, acceleration-controled translation model
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
        function obj = twoOrderModel(p_init, v_init, target)
            %TWOORDERMODEL constructor
            obj.p = p_init;
            obj.v = v_init;
            obj.target = target;

            % if recored the history
            if obj.isRecordHistory
                obj.history = zeros(twoOrderModel.setgetsteps(), 6);
                obj.recordOnce();
            end

            if isempty(twoOrderModel.setgetdt())
                % set the default dt
                twoOrderModel.setgetdt(0.1);
            end
            if isempty(twoOrderModel.setgetsteps())
                % set the default steps
                twoOrderModel.setgetsteps(2000);
            end
        end
        
        function obj = stepForward(obj, u)
            %stepForward update the state of the system, the dynamics of the model
            obj.v = obj.v + u * twoOrderModel.setgetdt();
            obj.p = obj.p + obj.v * twoOrderModel.setgetdt();

            if obj.isRecordHistory
                % record the state of the system
                obj.recordOnce();
            end
        end

        function obj = recordOnce(obj)
            %recordOnce record the state of the system
            %   recordOnce(obj) records the state of the system
            %   in the history array.
            obj.stepNow = obj.stepNow + 1;
            obj.history(obj.stepNow, :) = [obj.p', obj.v'];
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

    properties (SetAccess = private)
        % the safety control related parameters
        alpha1 = 0.1;
        alpha2 = 0.1;
        privArr = [];
    end

    methods
        function u_safe = distributedSafeFiltering(obj, u_ref, senseMat, privilege)
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
            [u_safe, modifValue, exitflag] = quadprog(H, f, CBFconsA, CBFconsb, [], [], [], [], [], obj.quadopt); % solve the quadratic programming problem

            % if the QP fails, use the minimum invasive result
            if exitflag < 0
                warning("QP error, using minimum invasive result")
                u_safe = twoOrderModel.minmaxlinear(CBFconsA, CBFconsb); 
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
                r_i = senseMat(i, 7);
                R_i = r_i + obj.r;
                
                d_i = norm(obj.p - p_i); % distance between the system and the sensed object
                n_hat = (obj.p - p_i) / d_i; % unitvector pointing from the sensed object to the system
                v_rel = obj.v - v_i; % relative velocity
                
                % control of agent i can guarentee safety by satisfying CBFf scaled by privilege
                % if n_hat' * v_rel >= 0
                %     v_priv = obj.privArr(i);
                % else
                %     v_priv = 1 - obj.privArr(i);
                % end
                CBFconsA(i, :) = -n_hat;
                CBFconsb(i) = (obj.alpha1 + obj.alpha2) * n_hat' * obj.v  ...
                                + obj.privArr(i)* ( obj.alpha1*obj.alpha2* (d_i - R_i) ...
                                    + 1/d_i * ( v_rel'*v_rel - n_hat'*v_rel * n_hat'*v_rel ) );

                if d_i < R_i - 1E-5
                    % give a warn of collision, collision may happen due to
                    % too large alpha and too large dt
                    warning("collision happened!");
                    d_i - R_i
                end
            end 
        end
    end


    %% optimization related
    properties
        quadopt = optimoptions("quadprog", "Display", "off");
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
        function plotHistory(obj, stepInterval)
            %plotHistory plot the history of the system with 3D quiver arrows
            %   plotHistory(obj, stepInterval) plots the trajectory of the system
            %   in 3D space with velocity vectors. stepInterval specifies the
            %   interval at which arrows are plotted (default is 30).
    
            if nargin < 2
                stepInterval = 40; % Default interval for plotting arrows
            end
    
            % Plot the trajectory
            p = plot3(obj.history(:, 1), obj.history(:, 2), obj.history(:, 3), 'LineWidth', 1.5, 'LineStyle', '--');
    
            % Prepare data for quiver arrows
            xs = obj.history(1:stepInterval:end, 1);
            ys = obj.history(1:stepInterval:end, 2);
            zs = obj.history(1:stepInterval:end, 3);
            vxs = obj.history(1:stepInterval:end, 4);
            vys = obj.history(1:stepInterval:end, 5);
            vzs = obj.history(1:stepInterval:end, 6);
    
            % Plot velocity vectors using quiver3
            % quiver3(xs, ys, zs, vxs, vys, vzs, 'AutoScale', 'on', 'Color', p.Color, 'LineWidth', 1.5, 'AutoScaleFactor', 1.2, 'MaxHeadSize', 20);
    
            % Highlight the target position
            plot3(obj.target(1), obj.target(2), obj.target(3), 'Marker', 'pentagram',  'MarkerSize', 10, 'LineWidth', 2, 'Color', p.Color);
   
        end
    end
end

