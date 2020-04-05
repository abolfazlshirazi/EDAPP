function [y,ceq] = G04_CON(x)
    % Problem:          G04
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    U = 85.334407+0.0056858*x(2)*x(5)+0.0006262*x(1)*x(4)-0.0022053*x(3)*x(5);
        y(1) = -U;
        y(2) = U-92;
    V = 80.51249+0.0071317*x(2)*x(5)+0.0029955*x(1)*x(2)+0.0021813*x(3)^2;
        y(3) = -V+90;
        y(4) = V-110;
    W = 9.300961+0.0047026*x(3)*x(5)+0.0012547*x(1)*x(3)+0.0019085*x(3)*x(4);
        y(5) = -W+20;
        y(6) = W-25;

    y=y'; 
    ceq = [];
end