function y = G07_OBJ(x)
    % Problem:          G07
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)

    y = x(1)^2+x(2)^2+x(1)*x(2)-14*x(1)-16*x(2)+(x(3)-10)^2+...
        4*(x(4)-5)^2+(x(5)-3)^2+2*(x(6)-1)^2+5*x(7)^2+...
        7*(x(8)-11)^2+2*(x(9)-10)^2+(x(10)-7)^2+45;
end