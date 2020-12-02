
function [strategy_success,strategy_prob]=strategy_prob_update(strategy_success,strategy_prob,num_strategy,K,n0,delta)

strategy_success(num_strategy)=strategy_success(num_strategy)+1;
for n=1:4
    strategy_prob(n)=( strategy_success(n)+n0)/(sum(strategy_success)+K*n0);  
end
if min(strategy_prob)<delta
    
    strategy_prob=ones(1,K)/K;
    strategy_success=zeros(1,K);
end

end