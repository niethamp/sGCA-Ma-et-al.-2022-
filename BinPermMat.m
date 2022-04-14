function [QueryTable,out] = BinPermMat(Groups)
%Groups=g
m=size(Groups,1);
QueryTable=[];

out=dec2bin(2^m-1:-1:0)-'0';

out=out(find(mean(out,2)~=0&mean(out,2)~=1),:); % remove invariant phenotypes

n_out=size(out,1);

QueryMat=zeros(n_out,max(Groups,[],'all'));

for i=1:n_out
    
    qname=cellstr(['P',num2str(i)]);
    
    for i2=1:m
        
        QueryMat(i,Groups(i2,1):Groups(i2,2))=out(i,i2);
 
        
    end
    
    QueryTable=[QueryTable;[array2table(qname),array2table(QueryMat(i,:))]];
    
    
end



end

