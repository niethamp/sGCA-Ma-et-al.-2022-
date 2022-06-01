function sGCA_Heatmap

%Centered colormap generator version 1.0.0 (1.85 KB) by Timothy Olsen
Default=colMapGen([0.87 0.11 0.76],[0.11 0.11 0.76],100);



%load Genetable with stats: Tab1 QueryMatrix, Tab2 Data with stats
[file,path]=uigetfile('','sGCA output','*.xlsx');
[~,name,~]=fileparts(file);
cd(path)
QueryTable=readtable([path,file],'Sheet','Query');
GeneTable=readtable([path,file],'Sheet','DataWithStats');
sourcedata=[path,file];


%User input
prompt = {'Maximal Correlation Distance', 'Minimal Fold-Regulation', 'Maximal Padj', 'Minimal MeanBase','Gene Annotation (-,*,+,GeneID, GeneSymbol; leave empty if none)','Labelsize','Print list (0/1)'};
dlgtitle = 'Thresholded sGCA heatmaps and gene lists (sGCA-Heatmap, v220601)';
dims = [1 100];
definput = {'0.25','2','0.05','20','','12','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

DistTh=str2num(answer{1});
FoldChangeTh=str2num(answer{2});
PvalTh=str2num(answer{3});
BaseTh=str2num(answer{4});
label=answer{5};
labelsize=str2num(answer{6});
printlist=str2num(answer{7});


Param=['[',answer{1},' ',answer{2},' ',answer{3},' ',answer{4},']'];

%datacol_start=


if isempty(label)==0
    %load Genetable with stats: Tab1 QueryMatrix, Tab2 Data with stats
    [file,path]=uigetfile('','Genes to be annotated','*.xlsx');
    LabledGenes=readtable([path,file]);
end

% filter table
GeneTableFilt=GeneTable(GeneTable.CorrelationDistance<DistTh&GeneTable.FoldChange>FoldChangeTh&GeneTable.Pval<PvalTh&GeneTable.MeanBase>BaseTh,:);

%sort table
for i=1:max(GeneTable.Phenotype)
    GeneTableFilt.PhenotypeCount(GeneTableFilt.Phenotype==i)=numel(GeneTableFilt.GeneID(GeneTableFilt.Phenotype==i));
end


sortGeneTableFilt=sortrows(GeneTableFilt,'MeanBase','descend'); %if there are no NNclusters
sortGeneTableFilt=sortrows(sortGeneTableFilt,'FoldChange','descend');
sortGeneTableFilt=sortrows(sortGeneTableFilt,'Pval','ascend');
sortGeneTableFilt=sortrows(sortGeneTableFilt,'CorrelationDistance','ascend');
sortGeneTableFilt=sortrows(sortGeneTableFilt,'Phenotype','ascend');
sortGeneTableFilt=sortrows(sortGeneTableFilt,'PhenotypeCount','descend');


% plot heatmap

sortGeneTableFiltNum=table2array(sortGeneTableFilt(:,vartype('numeric')));
sortGeneTableFiltNum=sortGeneTableFiltNum(:,1:(size(sortGeneTableFiltNum,2)-6)); %remove stat columns


figure,heatmap(sortGeneTableFiltNum,'Colormap',Default);
ax=gca;
ax.FontSize=labelsize;
ax.ColorScaling='scaledrows';
ylabels(1:size(sortGeneTableFiltNum,1))="";

if isempty(label)==0
    [IntersectList,IntersectInd]=intersect(sortGeneTableFilt.GeneID,LabledGenes.GeneID);



if label=="-"
    ylabels(IntersectInd)="-";
elseif label=="GeneID"
    ylabels(IntersectInd)=sortGeneTableFilt.GeneID(IntersectInd)+"-";
elseif label=="*"
    ylabels(IntersectInd)="*";
elseif label=="+"
    ylabels(IntersectInd)="+";
elseif label==">"
    ylabels(IntersectInd)=">";
elseif label=="GeneSymbol"
    ylabels(IntersectInd)=sortGeneTableFilt.GeneSymbol(IntersectInd)+"-";
end

end
ax=gca;
ax.YDisplayLabels=ylabels;
ax.XDisplayLabels=QueryTable.Properties.VariableNames(2:end);
grid off


% save data
writetable(sortGeneTableFilt,[[name,'_',Param,'_THtable'],'.xlsx']);
if printlist==1
writePhenotypeGenes3(GeneTable,[DistTh FoldChangeTh PvalTh BaseTh],[name,'_',Param,'_THlist']);
end
saveas(ax,[name,'_',Param,'_THheatmap.pdf']) %save heatmap as pdf
end




