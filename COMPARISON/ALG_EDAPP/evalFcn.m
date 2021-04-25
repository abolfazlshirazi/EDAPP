function varout=evalFcn(func,varin)
    n=size(varin,1); 
    varout=zeros(n,1);    
    
    for ii=1:n
        curvarout = func(varin(ii,:));
        varout(ii) = curvarout;
    end
end