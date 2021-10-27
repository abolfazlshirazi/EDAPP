function [x_output,FE_out] = extractSolution(Sol_arr,FE_arr,FEMax)
    for ii=1:numel(FE_arr)
        if FE_arr(ii) <= FEMax
            if ii==numel(FE_arr)
                disp('------ Info: ------');
                disp(['------ Total outer iterations: ' num2str(numel(FE_arr))])
                disp( '------ FE limit was not exceeded.' )
                disp('------------------');
                x_output = Sol_arr(ii,:);
                FE_out = FE_arr(ii);
                return;
            else
                continue;
            end
        else
            disp('------ Info: ------');
            disp(['------ Total outer iterations: ' num2str(numel(FE_arr))])
            disp(['------ FE limit exceeded at iteration: ' num2str(ii)])
            disp('------------------');
            x_output = Sol_arr(ii-1,:);
            FE_out = FE_arr(ii-1);
            return;
        end
    end
end
