% LLMexp.m
% the main script to cooperate LLM with distributer safety filter

clear all, clc, close all

dt_init = 0.5; % time interval of euler integration
N = 1000; % number of steps
n_agent = 6;

CWModel.setgetdt(dt_init);
CWModel.setgetsteps(N);

v_init = [0;0;0];

% starting position of agents
%      x    y    z
p1 =  [80,   0,    0;
       20,  -20,  -20;
       20,  20,   20;
       60,   20,   20;
       60,  -20,  -20;
       0,  0,    0];

% ending position of agents
%      x    y    z
p2 =  [0,   0,    0;
       60,   20,   20;
       60,  -20,  -20;
       20,  -20,  -20;
       20,  20,   20;
        80,  0,    0];

% add the LLM package
% install the package from
% https://github.com/matlab-deep-learning/llms-with-matlab.git
addpath("/home/chengrui/wk/llms-with-matlab/");

% the background prompt for the priority assignment task
backgroundPrompt = fileread("background/background.txt");

model = ollamaChat("gemma3:12b", ...
    "You are responsible for assigning priority for satellite swarms" + ...
     " this is the guideline you should follow: " + backgroundPrompt, ...
     "Temperature", 0.0);

% set the response format
formatStruct = struct("PriorityArray", ones(1, n_agent));

%% get priority assignment form LLM

% the question to LLM
form = 3;

switch form
    case 1
        question = "We now have 6 satellites, with satellite 1 mission critical. All other satellites are backup satellites. Give me the answer.";
    case 2
        question = "We now have 6 satellites, with satellite 1 mission critical. All other satellites are low fuel satellites. Give me the answer.";
    case 3
        question = "We now have 6 satellites, with all satellites are the same. Give me the answer.";
end

% generate the response
disp("getting result from LLM")
llmans = generate(model,question, "ResponseFormat", formatStruct);
llmans.PriorityArray

% get the priority matrix
priorityArr = llmans.PriorityArray;
priorityMatirx = zeros(n_agent, n_agent);
% construct priority matrix by weight of agent i over agent j
for i = 1:n_agent
    for j = i+1:n_agent
        priorityMatirx(i,j) = priorityArr(i) / (priorityArr(i) + priorityArr(j));
        priorityMatirx(j,i) = 1 - priorityMatirx(i,j);
    end
end

%% simulation
% init agent cell
agentCell = cell(1, n_agent);
for i = 1:n_agent
    agentCell{i} = CWModel(p1(i, :)', v_init, p2(i, :)');
end

firstArriveFlag = false;
% iterate for N steps, record the state of the system in history
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
        privilegeArray = priorityMatirx(j, :);
        privilegeArray(j) = [];
        u_safe = agentCell{j}.distributedSafeFiltering(u, senseMatAgentj, privilegeArray');
        U_safe((j-1)*3+1:(j-1)*3+3) = u_safe;
    end
    % step forward for each agent
    for j=1:length(agentCell)
        u_safe = U_safe((j-1)*3+1:(j-1)*3+3);
        agentCell{j}.stepForward(u_safe);
    end
    if norm(agentCell{1}.p - p2(1, :)') < 0.1 && ~firstArriveFlag
        disp("agent 1 arrived at step " + num2str(i))
        firstArriveFlag = true;
    end
end

% plot the trajectory of each agent
figure
hold on
for i = 1:length(agentCell)
    agentCell{i}.plotHistory();
end
hold off

view(3)
grid on
axis equal
xlim([0 80])
ylim([-30 30])
zlim([-30 30])