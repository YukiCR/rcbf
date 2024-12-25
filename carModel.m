classdef carModel < handle
    %% dynamics related
    properties
        % 小车的位置和方向
        R = 0.2; % 半径
        px = 0; % x 轴位置
        py = 0; % y 轴位置
        theta = 0; % 方向角
        target = [0,0,0]; % 目标点
        dt = 0.1; % 时间步长 (默认为 0.1)
        history = []; % 状态历史记录 (n x 3 矩阵)
        isRecordHistory = true;% 是否记录历史状态
        lastInput = []; % 上一个输入， 2x1
    end

    methods
        % 构造函数初始化位置、角度、时间步长和记录开关
        function obj = carModel(px_init, py_init, theta_init, target_init, dt_init, isRecordHistory_init)
            if nargin > 0
                obj.px = px_init;
                obj.py = py_init;
                obj.theta = theta_init;
                obj.target = target_init;
            end
            if nargin > 4
                obj.dt = dt_init;
            end
            if nargin > 5
                obj.isRecordHistory = isRecordHistory_init;
            end

            % 如果开启了历史记录功能，初始化时也记录状态
            if obj.isRecordHistory
                obj.recordState();
                obj.lastInput = [0;0];
            end
        end

        % 更新小车状态的方法
        function updateState(obj, in1, in2)
            if nargin == 3
                v = in1;
                omega = in2;
            end
            if nargin == 2 % vector form
                v = in1(1);
                omega = in1(2);
            end
            % v：速度，omega：角速度
            % 更新方向角
            obj.theta = obj.theta + omega * obj.dt;
            % 更新位置
            obj.px = obj.px + v * cos(obj.theta) * obj.dt;
            obj.py = obj.py + v * sin(obj.theta) * obj.dt;

            % 如果开启了历史记录功能，更新后记录状态
            if obj.isRecordHistory
                obj.recordState();
                obj.lastInput = [v; omega];
            end
        end

        % 获取当前状态的方法
        function [px, py, theta] = getState(obj)
            px = obj.px;
            py = obj.py;
            theta = obj.theta;
        end

        % 记录当前状态的方法
        function recordState(obj)
            % 将当前状态存入 history (一行数据：px, py, theta)
            state = [obj.px, obj.py, obj.theta];
            obj.history = [obj.history; state];
        end

        % 获取历史记录的方法
        function history = getHistory(obj)
            history = obj.history;
        end
    end

    %% nominal control related
    properties
        KP = [0.1, 1];
    end

    methods
        function [v,omega] = getPcontrol(obj)
            % see also 
            % CONTROL OF UNICYCLE TYPE ROBOTS
            % Tracking, Path Following and Point Stabilization
            errVec = [obj.target(1) - obj.px; obj.target(2) - obj.py];
            % 计算位置误差
            distance_error = errVec' * [cos(obj.theta); sin(obj.theta)];
            % 计算角度误差
            deltaTheta = wrapToPi(atan2(errVec(2), errVec(1)) - obj.theta);
            angle_error = deltaTheta;

            % 根据比例控制来计算线速度和角速度
            v = obj.KP(1) * distance_error;  % 线速度与位置误差成正比
            omega = obj.KP(2) * angle_error; % 角速度与角度误差成正比
        end
    end

    %% R-CBF related
    properties
        alpha = 0.5;
    end

    methods
        function filteredInput = RCBF_Filter(obj, agent, nominalInput)
            D = 0.1;
            % get states
            % position
            [pax, pay, patheta] = obj.getState();
            [pbx, pby, pbtheta] = agent.getState();
            Pa = [pax; pay];
            Pb = [pbx; pby];
            Pm = 0.5 * (Pa + Pb); % middle point
            % velocity
            % vOmegaA = obj.lastInput;
            % vOmegaB = agent.lastInput;
            % Va = vOmegaA(1) * [cos(patheta); sin(patheta)];
            % Vb = vOmegaB(1) * [cos(pbtheta); sin(pbtheta)];

            M = [cos(patheta), -D*sin(patheta);
                sin(patheta), D*cos(patheta)];

            na_hat = (Pa - Pm)/norm(Pa - Pm);
            B = 0.5*obj.alpha*(norm(Pa - Pb) - 2*obj.R);

            opt = optimoptions("quadprog","Display","off");

            % NOTE: THE BIG M REALLLLY MATTERS HERE
            % H = 2 * eye(2);
            % f = -2 * M * nominalInput; 
            % [filteredInput, ~, flag, ~] = quadprog(H,f,-na_hat',B,[],[],[],[],[],opt);
            % filteredInput = M \ filteredInput;

            H = 2 * M' * eye(2) * M;
            f = -2 * (M' * M) * nominalInput; 
            [filteredInput, ~, flag, ~] = quadprog(H,f,-na_hat' * M,B,[],[],[],[],[],opt);
            

            if flag < 0
                warning("QP error")
            end
        end
    end

    %% visualization related
    methods
        function plotTraj(obj)
            plot(obj.history(:,1), obj.history(:,2));
        end

        function plotHistory(obj)
            stepInterval = 30;  % 默认每stepInterval步绘制一次

            % 绘制轨迹
            p = plot(obj.history(1:stepInterval:end,1), obj.history(1:stepInterval:end,2), "LineWidth", 1, "LineStyle","--","Marker","o", "MarkerSize",8);

            xs = nan(length(1:stepInterval:length(obj.history)-1),1);
            ys = nan(length(1:stepInterval:length(obj.history)-1),1);
            vxs = nan(length(1:stepInterval:length(obj.history)-1),1);
            vys = nan(length(1:stepInterval:length(obj.history)-1),1);
            cnt = 1;

            for i = 1:stepInterval:length(obj.history)-1
                % 获取当前点的坐标和方向
                x = obj.history(i, 1);
                y = obj.history(i, 2);
                angle = obj.history(i, 3);

                % 计算速度矢量（假设速度为常数v）
                v = 1.0;  % 假设恒定速度为 1 m/s
                vx = v * cos(angle);  % x 方向速度
                vy = v * sin(angle);  % y 方向速度

                xs(cnt) = x;
                ys(cnt) = y;
                vxs(cnt) = vx;
                vys(cnt) = vy;

                cnt = cnt + 1;
            end

            % 使用quiver绘制速度箭头
            quiver(xs, ys, vxs, vys, 0.5, 'LineWidth', 1.0, 'MaxHeadSize', 1.0, "Color", p.Color);
        end

    end
end
