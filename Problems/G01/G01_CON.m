function [y,ceq] = G01_CON(x)
    % Problem:          G01
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y(1) = 2*x(1)+2*x(2)+x(10)+x(11)-10;
    y(2) = 2*x(1)+2*x(3)+x(10)+x(12)-10;
    y(3) = 2*x(2)+2*x(3)+x(11)+x(12)-10;
    y(4) = -8*x(1)+x(10);
    y(5) = -8*x(2)+x(11);
    y(6) = -8*x(3)+x(12);
    y(7) = -2*x(4)-x(5)+x(10);
    y(8) = -2*x(6)-x(7)+x(11);
    y(9) = -2*x(8)-x(9)+x(12);
    
    y=y';
    ceq = [];
end