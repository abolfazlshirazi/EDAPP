function [parent,smart] = evaluation(parent,smart,func)
    for ii=1:numel(parent)      
        parent{ii}.f_sam = evalFcn(func,parent{ii}.x_sam); 
    end
    for ii=1:numel(smart)
        smart{ii}.f_sam = evalFcn(func,smart{ii}.x_sam);
    end
end