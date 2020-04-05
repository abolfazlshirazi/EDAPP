function y = G09_OBJ(x)
    % Problem:          G09
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y = (x(1)-10)^2+5*(x(2)-12)^2+x(3)^4+3*(x(4)-11)^2+...
        10*x(5)^6+7*x(6)^2+x(7)^4-4*x(6)*x(7)-10*x(6)-8*x(7);
end