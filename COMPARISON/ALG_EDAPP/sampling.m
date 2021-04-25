function [parent,smart] = sampling(parent,smart,npop)
    percluster = floor(npop/(numel(parent)+numel(smart)));
    for ii=1:numel(parent)
        parent{ii}.x_sam = mvnrnd(parent{ii}.c,parent{ii}.sigma2,percluster);
    end
    for ii=1:numel(smart)
        smart{ii}.x_sam = mvnrnd(smart{ii}.c,smart{ii}.sigma2,percluster);
    end
end

