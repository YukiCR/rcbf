clear, clc, close all

%% init
formCell = cell(1, 4);

%% multi-sats, different forms

% set ingegral time interval and steps
dt_init = 0.5;
N = 800;
CWModel.setgetdt(dt_init);
CWModel.setgetsteps(N);
v_init = [0;0;0];

% generate agents on watch face
r = 30;
z1 = 10;
z2 = -10;
n_agent = 10;
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
% for i=2:n_agent
%     privilegeMat(1,i) = 0.9;
%     privilegeMat(i,1) = 0.1;
% end

% 1 for distributed with fixed P, 2 for centralized, 3 for distributed with optimized P,
% 4 for distributed assuming other agents are not controlled
form = 3;

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

                modif = modif + norm(u - u_safe);
            end
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
        end
end

% save results in formcell
formCell{form} = agentCell;


% =========== plotting =================
% plot the state of the system
labelFontSize = 12;
legendFontSize = 11;
width = 5;  
height = 5;  
xbound = [-30 30];
ybound = [-30 30];
zbound = [-25 25];

figure

% set
colororder("gem12");
C = colororder;

hold on 
for i = 1:length(agentCell)
    agentCell{i}.plotHistory('color', C(mod(i+1, length(C)), :));
end
hold off

% legend and labels
legend({'',  '$i$-th trajectory', '$i$-th starting direction', '$i$-th goal'}, 'Interpreter', 'latex','FontSize', legendFontSize);
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$y (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
zlabel('$z (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);

view(-25, 65)
axis equal
xlim(xbound)
ylim(ybound)
zlim(zbound)

% axes settigs
ax = gca;
ax.Box = "off";
ax.XLim = xbound;
ax.XGrid = "on";
ax.YGrid = "on";
ax.ZGrid = "on";
ax.FontSizeMode = 'manual';   % disable auto zoom
ax.FontSize = labelFontSize;
ax.FontName = "Times New Roman";
ax.XAxis.FontName = "Times New Roman";
ax.YAxis.FontName = "Times New Roman";
ax.ZAxis.FontName = "Times New Roman";

set(gca, 'LooseInset', get(gca,'TightInset'))
fig = gcf;
set(fig, 'Units', 'inches', 'Position', [1, 1, width, height]);

clear CWModel

%% compare distance of different cells
distanceArray = zeros(N, 4); % the minimum distance between agents in each form
for i = 1:4
    agentCell = formCell{i};
    r = agentCell{1}.r;
    for j = 1:N
        % the matrix of positions of agents, each row for an observation
        pMat = zeros(length(agentCell), 3); 
        for k = 1:length(agentCell)
            pMat(k, :) = agentCell{k}.history(j, 1:3);
        end
        distanceArray(j, i) = min(pdist(pMat)) - 2*r; % pdist computes the pairwise distance between rows of pMat
    end
end

% ===== plot the distance array =====
labelFontSize = 12;
legendFontSize = 11;
width = 6;  
height = 5;  
xbound = [-30 30];
ybound = [-30 30];
zbound = [-25 25];

% colororder
corder = [0.4660    0.6740    0.1880;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560];

figure
ax1 = axes();
ax2 = axes('Position', [0.5 0.30 0.35 0.2]);

% swicth to main axes
axes(ax1);
hold on
for i = 1:4
    plot((1:N)*dt_init, distanceArray(:, i), 'LineWidth', 2.0);
end
hold off

yline(0, "LineStyle","--", "LineWidth", 1.0)

axes(ax2);
slice = 250:400;
hold on
for i = 1:4
    plot(slice*dt_init, distanceArray(slice, i), 'LineWidth', 2.0);
end
hold off

yline(0, "LineStyle","--", "LineWidth", 1.0)

xlabel(ax1, 'Time (s)', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel(ax1, 'Minimum distance (m)', 'Interpreter', 'latex', 'FontSize', labelFontSize);
legend(ax1, {'Fixd priority',  'Centralized Control', 'Optimized priority', 'Non-cooperative'}, 'Interpreter', 'latex','FontSize', legendFontSize);

ax1.Box = "on";
ax1.XGrid = "on";
ax1.YGrid = "on";
ax1.YLim = [-1.5 19];
ax1.FontSizeMode = 'manual';   % disable auto zoom
ax1.FontSize = labelFontSize;
ax1.FontName = "Times New Roman";
ax1.XAxis.FontName = "Times New Roman";
ax1.YAxis.FontName = "Times New Roman";
ax1.ZAxis.FontName = "Times New Roman";

ax2.Box = "on";
ax2.XGrid = "on";
ax2.YGrid = "on";
ax2.GridLineStyle = "--";
ax2.LineWidth = 1.2;
ax2.YLim = [-1 1.5];
ax2.XLim = [slice(1) slice(end)]*dt_init;
ax2.FontSizeMode = 'manual';   % disable auto zoom
ax2.FontSize = labelFontSize;
ax2.FontName = "Times New Roman";
ax2.XAxis.FontName = "Times New Roman";
ax2.YAxis.FontName = "Times New Roman";
ax2.ZAxis.FontName = "Times New Roman";

colororder(corder)

set(gca, 'LooseInset', get(gca,'TightInset'))
fig = gcf;
set(fig, 'Units', 'inches', 'Position', [1, 1, width, height]);

%% distributed 2 agents

clear, clc

% set ingegral time interval and steps
dt_init = 0.2;
N = 1500;
CWModel.setgetdt(dt_init);
CWModel.setgetsteps(N);

p_init = [0;0;0];
v_init = [0;0;0];
target = [30;0;0];

agent1 = CWModel(p_init, v_init, target);
agent2 = CWModel(target, v_init, p_init + [0;0.01;0]);
agentCell = {agent1, agent2};

% set the privilege matrix
privform = 3;
switch privform
    case 1
        privilege = [0.5; 0.5];
    case 2
        privilege = [0.7; 0.3];
    case 3
        privilege = [1.0; 0.0];
end

% 1 for distributed, 2 for centralized
form = 1;

firstArriveFlag = false;
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

                modif = modif + norm(u - u_safe);
            end
            for j=1:length(agentCell)
                u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
                agentCell{j}.stepForward(u_safe);
            end
            if norm(agentCell{1}.p - target) < 0.1 && ~firstArriveFlag
                disp("agent 1 arrived at step " + num2str(i))
                firstArriveFlag = true;
            end
        end
end
% =========== plotting =================
% plot the state of the system
labelFontSize = 12;
legendFontSize = 11;
width = 5;  
height = 5;  
xbound = [-3 33];
ybound = [-11.5 5.5];
zbound = [-3.5 3.5];

figure

% set
colororder("gem12");
C = colororder;

hold on 
for i = 1:length(agentCell)
    agentCell{i}.plotHistory("stepInterval", 5, "quiverMultiplier", 5);
end
hold off

% legend and labels
legend({'',  '$i$-th trajectory', '$i$-th starting direction', '$i$-th goal'}, 'Interpreter', 'latex','FontSize', legendFontSize, 'Position', [0.60697354264693,0.297012949962734,0.39670016253436,0.123362747366192]);
xlabel('$x (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$y (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
zlabel('$z (\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', labelFontSize);

view(-0.1, 45)
axis equal
xlim(xbound)
ylim(ybound)
zlim(zbound)

% axes settigs
ax = gca;
ax.Box = "off";
ax.XLim = xbound;
ax.XGrid = "on";
ax.YGrid = "on";
ax.ZGrid = "on";
ax.FontSizeMode = 'manual';   % disable auto zoom
ax.FontSize = labelFontSize;
ax.FontName = "Times New Roman";
ax.XAxis.FontName = "Times New Roman";
ax.YAxis.FontName = "Times New Roman";
ax.ZAxis.FontName = "Times New Roman";

set(gca, 'LooseInset', [0 0 0 0])
fig = gcf;
set(fig, 'Units', 'inches', 'Position', [1, 1, width, height]);

clear CWModel