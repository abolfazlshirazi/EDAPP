function [x_sel,f_sel]=selection(x,f)
    trnc_factor = 0.5;
    npop = numel(f);
    nsel = round(trnc_factor*npop);

    if nsel < 2
        nsel = 2;
    end
    if nsel >= npop
        nsel = npop-1;
    end
        
    [f_sorted,ind_sorted] = sort(f);
    x_sorted = x(ind_sorted,:);
    
    x_sel = x_sorted(1:nsel,:);
    f_sel = f_sorted(1:nsel);
    
end