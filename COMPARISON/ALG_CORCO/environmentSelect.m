function  [p,objF,conV]=environmentSelect(p,objF,conV,trial,objFtrial,conVtrial,weights)

%
[popsize,n]=size(p);

% normalization
comfit=[objF;objFtrial];
comvoi=[conV;conVtrial];
normalcomfit=(comfit-min(comfit))/(max(comfit)-min(comfit)+1.e-15);
normalcomvoi=(comvoi-min(comvoi))/(max(comvoi)-min(comvoi)+1.e-15);

normalfit=normalcomfit(1:popsize); normalvoi=normalcomvoi(1:popsize);
normalfittrial=normalcomfit(1+popsize:end); normalvoitrial=normalcomvoi(1+popsize:end);

FIT=weights.*normalfit+(1-weights).*normalvoi;
FITtrial=weights.*normalfittrial+(1-weights).*normalvoitrial;

index=find(FITtrial < FIT);

p(index,:)=trial(index,:);
objF(index)=objFtrial(index);
conV(index)=conVtrial(index);