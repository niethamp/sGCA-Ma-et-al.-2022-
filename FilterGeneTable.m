function [filteredGeneTable, fGeneTableNum] = FilterGeneTable(GeneTable, GeneTableNum,ID_Col)


GeneTable.Properties.RowNames=GeneTable.GeneID;
GeneIDs=table2cell(GeneTable(:,ID_Col(1)));

[~,fGeneTableNum,fGeneIDs] = genelowvalfilter(GeneTableNum,GeneIDs);
[~,fGeneTableNum,fGeneIDs] = generangefilter(fGeneTableNum,fGeneIDs);
[~,fGeneTableNum,fGeneIDs] = genevarfilter(fGeneTableNum,fGeneIDs);

GeneIDs=array2table(fGeneIDs);
GeneIDs.Properties.VariableNames={'GeneID'};

filteredGeneTable=GeneTable(fGeneIDs,:);
end

