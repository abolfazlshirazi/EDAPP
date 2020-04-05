function y = G05_OBJ(x)
    % Problem:          G05
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y = 3*x(1)+1e-6*x(1)^3+2*x(2)+2e-6/3*x(2)^3;
end