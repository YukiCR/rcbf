clear, clc, close all

%%

% set ingegral time interval and steps
dt_init = 0.5;
N = 1500;
CWModel.setgetdt(dt_init);
CWModel.setgetsteps(N);

v_init = [0;0;0];

% generate agents on watch face
r = 30;
z1 = 10;
z2 = -10;
n_agent = 12;
theta = linspace(0, 2*pi, n_agent/2 + 1);
theta(end) = [];
p1 = [r*cos(theta); r*sin(theta); ones(1, n_agent/2) * z1];
p2 = [r*cos(theta + pi); r*sin(theta + pi); ones(1, n_agent/2) * z2];
agentCell = cell(1, n_agent);
Goals = zeros(3, n_agent);
for i = 1:n_agent/2
    agentCell{i} = CWModel(p1(:, i), v_init, p2(:, i));
    agentCell{i + n_agent/2} = CWModel(p2(:, i), v_init, p1(:, i));
    Goals(:, i) = p2(:, i);
    Goals(:, i + n_agent/2) = p1(:, i);
end

privilegeMat = (ones(n_agent) - eye(n_agent)) * 0.5;

% give large priority for agent 1
for i=2:n_agent
    privilegeMat(1,i) = 0.9;
    privilegeMat(i,1) = 0.1;
end

% 1 for distributed with fixed P, 2 for centralized, 3 for distributed with optimized P,
% 4 for distributed assuming other agents are not controlled
form = 1;

modif = 0;
switch form
    case 1
        % stepforward for N steps, record the state of the system in history
        for i = 1:N
            senseMat = zeros(length(agentCell), 7);
            for j = 1:length(agentCell)
                % get the state of the agents for sensing
                senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
            end
            U_safe = zeros(N*3, 1);
            for j = 1:length(agentCell)
                senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                u = agentCell{j}.getPIDcontrol();
                privilegeArray = privilegeMat(j, :);
                privilegeArray(j) = [];
                u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray');
                U_safe((j-1)*3+1:(j-1)*3+3) = u_safe;
                agentCell{j}.stepForward(u_safe);

                modif = modif + norm(u - u_safe);
            end
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end

    case 2
        % testing centralized form
        for i = 1:N
            senseMat = zeros(length(agentCell), 7);
            for j = 1:length(agentCell)
                % get the state of the agents for sensing
                senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
            end
            % get norminal control
            U_ref = zeros(length(agentCell)*3, 1);
            for j = 1:length(agentCell)
                u = agentCell{j}.getPIDcontrol();
                U_ref((j-1)*3+1:(j-1)*3+3) = u;
            end
            U_safe = CWModel.centralizedSafeFiltering(U_ref, senseMat);

            modif = modif + norm(U_ref - U_safe);
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end

    case 3
        % decentralized form with privilege optimization
        for i = 1:N
            senseMat = zeros(length(agentCell), 7);% randomized initialization
            for j = 1:length(agentCell)
                % get the state of the agents for sensing
                senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
            end
            U_safe = zeros(N*3, 1);
            for j = 1:length(agentCell)
                senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                u = agentCell{j}.getPIDcontrol();
                privilegeArray = privilegeMat(j, :);
                privilegeArray(j) = [];
                u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray');
                U_safe((j-1)*3+1:(j-1)*3+3) = u_safe;

                modif = modif + norm(u - u_safe);
            end
            % optimize every n steps
            if mod(i, 20) == 0
                % get the privilege matrix
                privilegeMat = CWModel.privilegeOptimize2(senseMat, privilegeMat, Goals)
            end
            % step forward
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end

    case 4
        % stepforward for N steps, record the state of the system in history
        for i = 1:N
            senseMat = zeros(length(agentCell), 7);
            for j = 1:length(agentCell)
                % get the state of the agents for sensing
                senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
            end
            U_safe = zeros(N*3, 1);
            for j = 1:length(agentCell)
                senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                u = agentCell{j}.getPIDcontrol();
                privilegeArray = privilegeMat(j, :);
                privilegeArray(j) = [];
                u_safe = agentCell{j}.distributedSafeFilteringBenchmark(u, senseMatAgentj, privilegeArray');
                U_safe((j-1)*3+1:(j-1)*3+3) = u_safe;
                agentCell{j}.stepForward(u_safe);

                modif = modif + norm(u - u_safe);
            end
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end
end

% plot the state of the system
figure

hold on 
for i = 1:length(agentCell)
    agentCell{i}.plotHistory();
end
hold off

view(3)
grid on
axis equal
xlim([-30 30])
ylim([-30 30])
zlim([-30 30])

clear CWModel

%% distributed 2 agents

clear, clc

% set ingegral time interval and steps
dt_init = 0.2;
N = 1500;
CWModel.setgetdt(dt_init);
CWModel.setgetsteps(N);

p_init = [0;0;0];
v_init = [0;0;0];
target = [50;0;0];

agent1 = CWModel(p_init, v_init, target);
agent2 = CWModel(target, v_init, p_init + [0;0.01;0]);
agentCell = {agent1};

% set the privilege matrix
privilege = [0.5; 0.5];


% 1 for distributed, 2 for centralized
form = 1;

modif = 0;
switch form
    case 1
        % stepforward for N steps, record the state of the system in history
        for i = 1:N
            senseMat = zeros(length(agentCell), 7);
            for j = 1:length(agentCell)
                % get the state of the agents for sensing
                senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
            end
            U_safe = zeros(N*3, 1);
            for j = 1:length(agentCell)
                senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
                u = agentCell{j}.getPIDcontrol();
                u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilege(j));
                U_safe((j-1)*3+1:(j-1)*3+3) = u_safe;
                agentCell{j}.stepForward(u_safe);

                modif = modif + norm(u - u_safe);
            end
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end
end

% plot the state of the system
figure

hold on 
for i = 1:length(agentCell)
    agentCell{i}.plotHistory();
end
hold off

view(3)
grid on
axis equal
xlim([0 50])
ylim([-20 20])
zlim([-20 20])

clear CWModel