function sGCA

% sGCA, v012922 by Philipp Niethammer: Interrogate mRNAseq count table for correlation with
% ideal phenotype profiles

version='v012922';

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
date=datestr(t);
name=[datestr8601(t),'sGCA'];
warning off

AllOutput=[];
%User input
prompt = {'GeneID columns','Filter counts (0=no,1=yes)' ,'RawCount scaling (median of ratios; 0=no, 1=yes)','Groups:'};
dlgtitle = 'simple Gene Correlation Analysis by P.Niethammer (sGCA, v012922)';
dims = [1 100];
definput = {'1,2','1','0','1,3;4,6;7,9;10,12;13,15;16,18'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


ID_Col=str2num(answer{1});
Filt=str2num(answer{2});
CountNorm=str2num(answer{3});
g=str2num(answer{4});
a=1;


% Import mRNAseq count table from Excel
[file,path]=uigetfile('','mRNAseq counts','*.xlsx');
cd(path)
GeneTable=readtable([path,file]);
GeneTable=GeneTable(:,min(ID_Col):size(GeneTable,2));
sourcedata=[path,file]

settings.sourcedata=sourcedata;
settings.date=date;
settings.version=version;

% prune table to selected groups
ind_start=1;
g_ind_new=[];
ind_end=0;
g_ind=[];
g_new=[];
for n_group=1:size(g,1)
    ind=g(n_group,1):g(n_group,2);
    ind_end=(ind_start+(max(ind)-min(ind)));
    ind_new=ind_start:ind_end;
    g_ind=[g_ind, ind]; 
    g_ind_new=[g_ind_new,ind_new];
    g_new=[g_new;[ind_start,ind_end]];
    ind_start=ind_end+1;
end

GeneTable=GeneTable(:,[ID_Col,(g_ind+max(ID_Col))]);

if CountNorm==1
    GeneTableNum=NormalizeCounts(table2array(GeneTable(:,vartype('numeric'))));
else
    GeneTableNum=table2array(GeneTable(:,vartype('numeric')));
end



GeneTable(:,max(ID_Col)+1:size(GeneTable,2))=array2table(GeneTableNum);
NumberInputGenes=size(GeneTable,1)

if Filt==1
[GeneTable, GeneTableNum]=FilterGeneTable(GeneTable, GeneTableNum,ID_Col);
end

GenesRemovedByFilter=NumberInputGenes-size(GeneTable,1)


    % Generate permutation query table
    [QueryTable,GroupPermut]=BinPermMat(g_new);
    QueryTable.Properties.VariableNames(1)={'Phenotype'};
    QueryTable.Properties.VariableNames(2:size(QueryTable,2))=GeneTable.Properties.VariableNames(numel(ID_Col)+1:size(GeneTable,2));


QueryTableNum=table2array(QueryTable(:,vartype('numeric')));



% Calculate correlation distances between gene expression profiles & query profiles 

PhenoDist=pdist2(GeneTableNum, QueryTableNum,'correlation');
  

for NoP=1:size(QueryTableNum,1)
    NoP
  Treated=GeneTableNum(:,find(QueryTableNum(NoP,:)==1));
  Untreated=GeneTableNum(:,find(QueryTableNum(NoP,:)==0),:);
  
  meanTreated=mean(Treated,2 );
  meanUntreated=mean(Untreated,2);

  % compute the mean and the log2FC
    meanBase(:,NoP) = (meanTreated + meanUntreated) / 2;
    foldChange(:,NoP) = meanTreated ./ meanUntreated;
    log2FC(:,NoP)= log2(foldChange(:,NoP));
    tLocal = nbintest(Treated,Untreated,'VarianceLink','LocalRegression');
    figure(1),plotVarianceLink(tLocal)
    padj(:,NoP)= mafdr(tLocal.pValue,'BHFDR',true);
end

%Rank phenotypes by correlation distance

[~,Phenotype]=min(PhenoDist');
CorrelationDistance=zeros(size(GeneTableNum,1),1);
FoldChange=zeros(size(GeneTableNum,1),1);
Pval=zeros(size(GeneTableNum,1),1);
MeanBase=zeros(size(GeneTableNum,1),1);

for i=1:size(GeneTableNum,1)
CorrelationDistance(i)=PhenoDist(i,Phenotype(i));
MeanBase(i)=meanBase(i,Phenotype(i))';
FoldChange(i)=foldChange(i,Phenotype(i));
Pval(i)=padj(i,Phenotype(i));
end

CorrelationDistance=array2table(CorrelationDistance);
FoldChange=array2table(FoldChange);
Pval=array2table(Pval);
MeanBase=array2table(MeanBase);
Phenotype=Phenotype';


% save data
Phenotype=array2table(Phenotype);
writetable(QueryTable,[name,'.xlsx'],'sheet','Query');
exGeneTable=[GeneTable,Phenotype,CorrelationDistance, FoldChange,Pval,MeanBase];
writetable(exGeneTable,[name,'.xlsx'],'sheet','DataWithStats','WriteMode','append');
save(name);

close all
end