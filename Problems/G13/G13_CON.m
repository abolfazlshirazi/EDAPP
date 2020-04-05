function [y,ceq] = G13_CON(x)
    % Problem:          G13
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = abs(sum(x.^2)-10)-1e-3;
    y(2) = abs(x(2)*x(3)-5*x(4)*x(5))-1e-3;
    y(3) = abs(x(1)^3+x(2)^3+1)-1e-3;

    y=y'; 
    ceq = [];
end