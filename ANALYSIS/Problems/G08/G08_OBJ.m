function y = G08_OBJ(x)
    % Problem:          G08
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y = -(sin(2*pi*x(1))^3*sin(2*pi*x(2)))/(x(1)^3*(x(1)+x(2)));
end