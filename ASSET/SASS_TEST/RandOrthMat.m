function R = RandOrthMat(n, t)
    % orthogonal matrix approx
    R = eye(n);
    l = randperm(n);
    % t = 1e-8;
    for ii = 1:floor(n/2)
        i = 2*(ii-1)+1;
        R(l(i),l(i)) = sin(t);
        R(l(i+1),l(i+1)) = sin(t);
        R(l(i),l(i+1)) = cos(t);
        R(l(i+1),l(i)) = -cos(t);
    end 
end 