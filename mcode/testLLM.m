clear all, clc, close all

% add the LLM package
addpath("/home/chengrui/wk/llms-with-matlab/");

backgroundPrompt = fileread("background/background.txt");
model = ollamaChat("gemma3:12b", ...
    "You are responsible for assigning priority for satellite swarms" + ...
     " this is the guideline you should follow: " + backgroundPrompt, ...
     "Temperature", 0.0);
question = "We now have 6 satellites, with satellite 1 mission critical. satellite 3 low fuel. And satellite 6 is backup satellite. Give me the answer.";
formatStruct = struct("PriorityArray", ones(1, 6));
generate(model,question, "ResponseFormat", formatStruct)
