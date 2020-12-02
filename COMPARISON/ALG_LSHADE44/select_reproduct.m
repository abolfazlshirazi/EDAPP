function y = select_reproduct(indiv,s,lu,D,CR,F,population,popsize,num_strategy)
%% DE operator %%
y = de_crossover(indiv, s,population, popsize,D, CR, F ,num_strategy);
%% repair operator %%
y = min(max(y,lu(1,:)),lu(2,:));
end

function y=de_crossover(indiv,s, population,popsize, D, CR, F,num_strategy)
%% parent selection %%
rand_index = randi(popsize,1,2);
while rand_index(1) == rand_index(2)||rand_index(1) == s||rand_index(2) == s
    rand_index = randi(popsize,1,2);
end
parent1 = population(rand_index(1),:);
parent2 = population(rand_index(2),:);

if num_strategy<3
p=(0.2-1/popsize)*rand+1/popsize; %get the pbest;
pnum = floor(p*popsize);
population=epsilon_sort(population,0);
pbest = population(randi(pnum),:);
indiv=indiv(1:D);
parent1=parent1(1:D);
parent2=parent2(1:D);
pbest=pbest(1:D);
end




    
%% strategy 1: DE/current to pbest/1/bin
if num_strategy==1
%%%%% current to pbest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = indiv+ F * (pbest-indiv) + F * (parent1 - parent2);
%%%%%% binomial crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cross_rand = rand(1,D);
    cross_rand(randi(D)) = 0; % at least one decision component to crossover
    cross_ID=cross_rand < CR;
    y(~cross_ID)= indiv(~cross_ID);
end

%% strategy 2: DE/current to pbest/1/exp
if num_strategy==2
%%%%% current to pbest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = indiv+ F * (pbest-indiv) + F * (parent1 - parent2);
%%%%%% exponent crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cross_ID=false(1,D);
    start_ID=randi(D);
    place_ID=start_ID;
    cross_ID(place_ID)=true;
    while rand<CR^(place_ID-1) && sum(cross_ID)<D
        place_ID=place_ID+1;
        if place_ID <= D
            cross_ID(place_ID)=true;
        else
            place_ID=mod(place_ID,D);
            cross_ID(place_ID)=true;
        end
    end
    y(~cross_ID)= indiv(~cross_ID);
end

%% strategy 3: DE/randr1/1/bin
if num_strategy==3
    %%%%%% replace r1, r2 and r3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r1=indiv;
    r2=parent1;
    r3=parent2;
    if (r1(end)>r2(end))||((r1(end)==r2(end))&&(r1(end-1)>r2(end-1)))
        temp=r1;
        r2=temp;
        r1=r2;
    end
    if (r1(end)>r3(end))||((r1(end)==r3(end))&&(r1(end-1)>r3(end-1)))
        temp=r1;
        r3=temp;
        r1=r3;
    end
    r1=r1(1:D);
    r2=r2(1:D);
    r3=r3(1:D);
    y=r1+F*(r2-r3);
%%%%%% binomial crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cross_rand = rand(1,D);
    cross_rand(randi(D)) = 0; % at least one decision component to crossover
    cross_ID=cross_rand < CR;
    y(~cross_ID)= r1(~cross_ID);
end


%% strategy 4: DE/randr1/1/exp
if num_strategy==4
%%%%%% replace r1, r2 and r3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r1=indiv;
    r2=parent1;
    r3=parent2;
    if (r1(end)>r2(end))||((r1(end)==r2(end))&&(r1(end-1)>r2(end-1)))
        temp=r1;
        r2=temp;
        r1=r2;
    end
    if (r1(end)>r3(end))||((r1(end)==r3(end))&&(r1(end-1)>r3(end-1)))
        temp=r1;
        r3=temp;
        r1=r3;
    end
    r1=r1(1:D);
    r2=r2(1:D);
    r3=r3(1:D);
    y=r1+F*(r2-r3);
%%%%%% exponent crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cross_ID=false(1,D);
    start_ID=randi(D);
    place_ID=start_ID;
    cross_ID(place_ID)=true;
    while rand<CR^(place_ID-1) && sum(cross_ID)<D
        place_ID=place_ID+1;
        if place_ID <= D
            cross_ID(place_ID)=true;
        else
            place_ID=mod(place_ID,D);
            cross_ID(place_ID)=true;
        end
    end
    y(~cross_ID)= r1(~cross_ID);
end



end