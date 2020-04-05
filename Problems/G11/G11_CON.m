function [y,ceq] = G11_CON(x)
    % Problem:          G11
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = abs(x(2)-x(1)^2)-1e-3;
    
    y=y'; 
    ceq = [];
end