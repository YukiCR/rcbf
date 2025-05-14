clear, clc, close all

% TODO: centralized realize

% set ingegral time interval and steps
dt_init = 0.5;
N = 2000;
twoOrderModel.setgetdt(dt_init);
twoOrderModel.setgetsteps(N);

p_init = [0;0;0];
v_init = [0;0;0];
target = [50; 0; 0];

agent1 = twoOrderModel(p_init, v_init, target);
agent2 = twoOrderModel(target, v_init, p_init + [0; 0.01; 0]);
agent3 = twoOrderModel([40; 10; 10], v_init, [0; -10; -10]);
agent4 = twoOrderModel([0; -10; -10], v_init, [40; 10; 10]);

agentCell = {agent1, agent2, agent3, agent4};
privilege = [0.4; 0.4; 0.4; 0.5];

privilegeMat = [0,   0.4, 0.4, 0.5;
                0.6, 0,   0.5, 0.3;
                0.6, 0.5, 0,   0.2;
                0.5, 0.7, 0.8, 0.0];

% stepforward for N steps, record the state of the system in history
for i = 1:N
    senseMat = zeros(length(agentCell), 7);
    for j = 1:length(agentCell)
        % get the state of the agents for sensing
        senseMat(j, :) = [agentCell{j}.p', agentCell{j}.v', agentCell{j}.r];
    end
    for j = 1:length(agentCell)
        senseMatAgentj = senseMat((1:length(agentCell)~=j), :);
        u = agentCell{j}.getPIDcontrol();
        privilegeArray = privilegeMat(:, j);
        privilegeArray(j) = [];
        u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray);
        agentCell{j}.stepForward(u_safe);
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

clear twoOrderModel