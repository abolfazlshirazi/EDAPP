function CR=CR_product(CR0)
CR=randn*0.1+CR0;

while CR>1||CR<=0
    CR=randn*0.1+CR0;
end

end