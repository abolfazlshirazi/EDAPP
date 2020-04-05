function y = G02_OBJ(x)
    % Problem:          G02
    % Function Type:    Objectives
    % Date:             Nov. 2019
    % By:               Abolfazl Shirazi (ashirazi@bcamath.org)
    
    N = 20; 
    totalsum = 0;
    for ii=1:N
        totalsum = totalsum+ii*x(ii)^2;
    end
    y = -abs((sum(cos(x).^4)-2*prod(cos(x).^2))/sqrt(totalsum));
end