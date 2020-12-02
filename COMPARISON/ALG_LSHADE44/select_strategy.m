function num_strategy=select_strategy(strategy_prob)

%% select strategy %%
rand_strategy=rand(1);
if rand_strategy<strategy_prob(1)
    num_strategy=1;
elseif rand_strategy<strategy_prob(1)+strategy_prob(2)
    num_strategy=2;
elseif rand_strategy<strategy_prob(1)+strategy_prob(2)+strategy_prob(3)
    num_strategy=3;
else
    num_strategy=4;
end
end