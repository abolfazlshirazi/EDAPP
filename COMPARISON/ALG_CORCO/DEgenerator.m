function trial=DEgenerator(p,objF,conV,archive,weights,minVar,maxVar)

% A. Shirazi
% global gen maxGen problemSetNum
problemSetNum = 2006; % arbitrary!

global gen maxGen
lu=[minVar;maxVar];
[popsize,n]=size(p);

trial=zeros(popsize,n);

normalvoi=(conV-min(conV))./(max(conV)-min(conV)+1.e-15);
normalfit=(objF-min(objF))./(max(objF)-min(objF)+1.e-15);

%%
for i=1:popsize

    l=rand;
    if l < 1/2
        F  = .6;
    elseif l < 2/3
        F  = .8;
    else
        F = 1;
    end
    
     l=rand;
    if l < 1/2
        CR  = .1;
    elseif l < 2/3
        CR  = .2;
    else
        CR = 1.0;
    end

    %% Initialization for CR
    indexset=1:popsize;
    indexset(i)=[];
    r1=floor(rand*(popsize-1))+1;
    xr1=indexset(r1);
    indexset(r1)=[];
    r2=floor(rand*(popsize-2))+1;
    xr2=indexset(r2);
    indexset(r2)=[];
    
    r=floor(rand*(popsize-3))+1;
    target=p(indexset(r),:);

    if rand < 0.5
        FIT=weights(i)*normalfit+(1-weights(i))*normalvoi;
        %FIT=weights*normalfit+(1-weights)*normalvoi;
        %FIT=normalfit;
        [~,best]=min(FIT);
        v=p(xr1,:)+F*(p(best,:)-p(xr1,:))+F*(target-p(xr2,:));
        flag=0;
        
    else
        r3=floor(rand*(popsize-3))+1;
        xr3=indexset(r3);
        v=p(i,:)+rand*(p(xr1,:)-p(i,:))+F*(p(xr3,:)-p(xr2,:));
        %v=p(i,:)+F*(p(xr3,:)-p(xr2,:));
        flag=1; 
        
    end
       
    w = find(v < lu(1, :));
    if ~isempty(w)
        l=rand;
        if l <= 0.5
            v(1, w) = 2 * lu(1, w) -  v(1, w);
            w1 = find( v(1, w) > lu(2, w));
            if ~isempty(w1)
                v(1, w(w1)) = lu(1, w(w1));
            end
        else
            if gen<0.5*maxGen && problemSetNum == 2006
               v(1, w) =  lu(1, w);
            else
               v(1, w) =  lu(2, w);
            end
        end
    end
    
    y = find(v > lu(2, :));
    if ~isempty(y)
        l=rand;
        if l <=0.7
            v(1, y) =  2 * lu(2, y) - v(1, y);
            y1 = find(v(1, y) < lu(1, y));
            if ~isempty(y1)
                v(1, y(y1)) = lu(2, y(y1));
            end
         else
             
            if gen<0.5*maxGen && problemSetNum == 2006
                v(1, y) =  lu(2, y);
            else
                v(1, y) =  lu(1, y);
            end
        end
    end
    
    if flag==0
        % Binomial crossover
        t = rand(1, n) < CR;
        j_rand = floor(rand * n) + 1;
        t(1, j_rand) = 1;
        t_ = 1 - t;
        trial(i, :) = t .* v + t_ .* p(i, :);
    else
        trial(i,:)=v;
    end
end