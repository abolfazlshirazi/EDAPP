function [x_new,f_new,x_best,f_best] = replacement(parentclusters,smartclusters,x,f)
    % Best Elitism    
    npop = numel(f);

        x_sam = [];
        f_sam = [];
        for jj=1:numel(parentclusters)
            x_sam = [x_sam;parentclusters{jj}.x_sam];
            f_sam = [f_sam;parentclusters{jj}.f_sam];
        end
        for jj=1:numel(smartclusters)
            x_sam = [x_sam;smartclusters{jj}.x_sam];
            f_sam = [f_sam;smartclusters{jj}.f_sam];
        end
        
    x = [x;x_sam];
    f = [f;f_sam];
    
    [f_sorted,ind]= sort(f);
    x_sorted = x(ind,:);
    
    x_new = x_sorted(1:npop,:);
    f_new = f_sorted(1:npop,:);   
    
    x_best = x_new(1,:);
    f_best = f_new(1);
end