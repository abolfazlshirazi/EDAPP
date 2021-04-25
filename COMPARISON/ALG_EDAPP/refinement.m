function [parentclusters,num_of_ref] = refinement(parentclusters,x_sel,f_sel)
    num_of_ref = 0;
    for ii=1:numel(parentclusters)
        if size(parentclusters{ii}.x,1)==1    
            num_of_ref = num_of_ref+1;

            tocheck = (x_sel == repmat(parentclusters{ii}.x,size(x_sel,1),1));

            thisind = 0;
            for kk=1:size(tocheck,1)
                if tocheck(kk,:) == ones(1,size(x_sel,2))
                    thisind = kk;
                    break;
                end
            end
            if thisind~=0
                x_rem = x_sel((1:size(x_sel,1))'~=thisind,:);
                f_rem = f_sel((1:size(x_sel,1))'~=thisind);
            else

                x_rem = x_sel;
                f_rem = f_sel;
            end

            d = x_rem - repmat(parentclusters{ii}.x,size(x_rem,1),1);
            [~,q2] = min(sum(d.^2,2));
            thepointtoaddtothiscluster = x_rem(q2,:);
            thepointtoaddtothisclusterf = f_rem(q2);

            parentclusters{ii}.x = [parentclusters{ii}.x; thepointtoaddtothiscluster];
            parentclusters{ii}.f = [parentclusters{ii}.f; thepointtoaddtothisclusterf];
        end
    end
end