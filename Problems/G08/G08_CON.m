function [y,ceq] = G08_CON(x)
    % Problem:          G08
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = x(1)^2-x(2)+1;
    y(2) = 1-x(1)+(x(2)-4)^2;

    y=y'; 
    ceq = [];
end