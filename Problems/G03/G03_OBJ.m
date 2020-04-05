function y = G03_OBJ(x)
    % Problem:          G03
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    N = length(x);
    y = -sqrt(N)^N*prod(x);
end