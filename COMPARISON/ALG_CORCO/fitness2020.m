function [objF, conV]=fitness2020(x, prob_k)

global const_num;

[f,g,h] = cec20_func(x,prob_k);

% Obtain the fitness
objF= f;
conV = (sum([sum(g.*(g>0),1); sum(abs(h).*(abs(h)>1e-4),1)])./const_num);
conV = conV';
