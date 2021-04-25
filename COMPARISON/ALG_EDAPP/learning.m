function [parent,smart] = learning(parent,smart)
    for ii=1:numel(parent)
        parent{ii}.sigma2 = var(parent{ii}.x);
    end
    for ii=1:numel(smart)
        [~,bestind] = min(smart{ii}.f);
        
        smart{ii}.c = smart{ii}.x(bestind,:); 
        
        smart{ii}.sigma2 = 2.*abs(smart{ii}.c - smart{ii}.parent);
    end
end
