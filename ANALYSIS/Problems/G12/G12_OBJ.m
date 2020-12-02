function y = G12_OBJ(x)
    % Problem:          G12
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    y = -(1-0.01*((x(1)-5)^2+(x(2)-5)^2+(x(3)-5)^2));
end