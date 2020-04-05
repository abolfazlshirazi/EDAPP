function y = G01_OBJ(x)
    % Problem:          G01
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)

    x1 = x(1:4);
    x2 = x(5:13);
    y = 5*sum(x1)-5*sum(x1.*x1)-sum(x2);
end