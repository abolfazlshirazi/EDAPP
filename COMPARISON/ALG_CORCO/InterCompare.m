function [size1,size2]=InterCompare(conleadpop_obj,conleadpop_con,objleadpop_obj,objleadpop_con)


mix_obj = [objleadpop_obj;conleadpop_obj];
mix_con = [conleadpop_con;objleadpop_con];

[~,objIndex]=sort(mix_obj);
[~,conIndex]=sort(mix_con);

popsize = size(conleadpop_obj,1);
temp1 = find(objIndex(1:popsize)>popsize);
temp2 = find(conIndex(1:popsize)>popsize);
size1 = size(temp1,1);
size2 = size(temp2,1);