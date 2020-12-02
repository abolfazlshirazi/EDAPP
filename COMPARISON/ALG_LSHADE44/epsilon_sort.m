function population=epsilon_sort(population,epsilon)
    temp_population=population;
    temp_population(temp_population(:,end)<=epsilon,end)=0;
    [~,rank]=sort(temp_population(:,end));
    fea=temp_population(population(:,end)==0);
    feasize=size(fea,1);
    if feasize~=0
    fea_rank = rank(1:feasize);
    fea_obj=temp_population(fea_rank,end-1);
    [~,fea_rank2]=sort(fea_obj);
    fea_rank=fea_rank(fea_rank2);
    final_rank=[fea_rank;rank(feasize+1:end)];
    else
    final_rank=rank;
    end
    population=population(final_rank,:);
end