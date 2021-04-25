function [parent,smart] = repairing(parent,smart,lb,ub)
    for ii=1:numel(parent)
        n = size(parent{ii}.x_sam,1);
        parent{ii}.x_sam = max(repmat(lb,n,1),parent{ii}.x_sam);
        parent{ii}.x_sam = min(repmat(ub,n,1),parent{ii}.x_sam);        
    end
    for ii=1:numel(smart)
        n = size(smart{ii}.x_sam,1);
        smart{ii}.x_sam = max(repmat(lb,n,1),smart{ii}.x_sam);
        smart{ii}.x_sam = min(repmat(ub,n,1),smart{ii}.x_sam);        
    end
end
