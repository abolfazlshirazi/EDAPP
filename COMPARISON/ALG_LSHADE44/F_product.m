function F=F_product(F0)
F=cauchy_rand(F0,0.1);
while F<=0 ||F>=1
    F=cauchy_rand(F0,0.1);
end

end