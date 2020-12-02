function [y,ceq] = G10_CON(x)
    % Problem:          G10
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)

    y(1) = -1+0.0025*(x(4)+x(6));
    y(2) = -1+0.0025*(-x(4)+x(5)+x(7));
    y(3) = -1+0.01*(-x(5)+x(8));
    y(4) = 100*x(1)-x(1)*x(6)+833.33252*x(4)-83333.333;
    y(5) = x(2)*x(4)-x(2)*x(7)-1250*x(4)+1250*x(5);
    y(6) = x(3)*x(5)-x(3)*x(8)-2500*x(5)+1250000;

    y=y';
    ceq = [];
end