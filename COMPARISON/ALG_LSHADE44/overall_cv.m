function result = overall_cv(g,h)
h=abs(h)-1e-4;
cv=[g,h];
cv(cv < 0) = 0;
result = sum(cv,2);
end


