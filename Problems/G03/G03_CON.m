function [y,ceq] = G03_CON(x)
    % Problem:          G03
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = abs(sum(x.^2)-1)-1e-3;

    y=y'; 
    ceq = [];
end