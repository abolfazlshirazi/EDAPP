function [y,ceq] = G05_CON(x)
    % Problem:          G05
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = x(3)-x(4)-0.55;
    y(2) = x(4)-x(3)-0.55;
    y(3) = abs(1000*(sin(-x(3)-0.25)+sin(-x(4)-0.25))+894.8-x(1))-1e-3;
    y(4) = abs(1000*(sin(x(3)-0.25)+sin(x(3)-x(4)-0.25))+894.8-x(2))-1e-3;
    y(5) = abs(1000*(sin(x(4)-0.25)+sin(x(4)-x(3)-0.25))+1294.8)-1e-3;

    y=y';
    ceq = [];
end

