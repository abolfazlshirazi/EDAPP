function [y,cq] = G02_CON(x)
    % Problem:          G02
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    n = 20;
    y(1) = -prod(x)+0.75;
    y(2) = sum(x)-7.5*n;
    
    y=y';
    cq = [];
end