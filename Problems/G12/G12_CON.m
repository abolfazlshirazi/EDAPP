function [y,ceq] = G12_CON(x)
    % Problem:          G12
    % Function Type:    Constraints
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    for P=1:9
        for Q=1:9
            for R=1:9
                z(P,Q,R) = (x(1)-P)^2+(x(2)-Q)^2+(x(3)-R)^2-0.0625;
            end
        end
    end
    for P=1:9
        for Q=1:9
            Z1(P,Q) = min(z(P,Q,:));    
        end
        Z2(P) = min(Z1(P,:));
    end
    y = min(Z2);
    ceq = [];
end