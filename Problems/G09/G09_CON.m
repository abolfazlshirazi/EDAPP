function [y,ceq] = G09_CON(x)
    % Problem:          G09
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)

    V1 = 2*x(1)^2;
    V2 = x(2)^2;
    y(1) = V1+3*V2^2+x(3)+4*x(4)^2+5*x(5)-127;
    y(2) = 7*x(1)+3*x(2)+10*x(3)^2+x(4)-x(5)-282;
    y(3) = 23*x(1)+V2+6*x(6)^2-8*x(7)-196;
    y(4) = 2*V1+V2-3*x(1)*x(2)+2*x(3)^2+5*x(6)-11*x(7);

    y=y';
    ceq = [];
end