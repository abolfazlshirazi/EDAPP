function [y,ceq] = G06_CON(x)
    % Problem:          G06
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = -(x(1)-5)^2-(x(2)-5)^2+100;
    y(2) = (x(1)-6)^2+(x(2)-5)^2-82.81;

    y=y';
    ceq = [];
end