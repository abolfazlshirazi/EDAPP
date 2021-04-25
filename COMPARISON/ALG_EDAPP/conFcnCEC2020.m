function out = conFcnCEC2020(x)
    global fcnNum const_num;
    [~,g,h] = cec20_func(x,fcnNum);    
    out = (sum([sum(g.*(g>0),1); sum(abs(h).*(abs(h)>1e-4),1)])./const_num);
end