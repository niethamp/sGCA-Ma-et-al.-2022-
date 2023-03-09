function sGEA_230224

% note: this is a beta version. Errors are possible. If you find one,
% please notify me niethamp@mskcc.org

version='v230224';

%Centered colormap generator version 1.0.0 (1.85 KB) by Timothy Olsen, not

Default=colMapGen([0.87 0.11 0.76],[0.11 0.11 0.76],100);

warning off

%User input
prompt = {'Analyzed Genes'};
dlgtitle = 'simple Gene Enrichment Analysis by P.Niethammer (sGEA, v230224)';
dims = [1 100];
definput = {'25000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
nAllGenes=str2num(answer{1});

%load marker gene list
[fileMarker,path]=uigetfile('','Marker genes (e.g., filtered sGCA-table)','*.xlsx');
[~,name,~]=fileparts(fileMarker);
cd(path)
sourcedataMarker=[path,fileMarker];
MarkerTable=readtable(sourcedataMarker);

%load GenetableGCA
[fileGCA,path]=uigetfile('','Filtered sGCA-table','*.xlsx');
[~,name,~]=fileparts(fileGCA);
cd(path)
sourcedataGCA=[path,fileGCA];
GeneTableGCA=readtable(sourcedataGCA);

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
date=datestr(t);
name=[datestr8601(t),'sGEA'];

fileResult=[fileGCA,'_vs_',fileMarker];

GcaClusters=unique(GeneTableGCA.Phenotype,'stable');
MarkerClusters=unique(MarkerTable.Phenotype, 'stable');


for nIP=1:numel(GcaClusters)

    GcaClusterGenes=GeneTableGCA.GeneSymbol(GeneTableGCA.Phenotype==GcaClusters(nIP));
    
    RowName(nIP)=strcat("IP",num2str(GcaClusters(nIP)));

    for nIP2=1:numel(MarkerClusters)

        if nIP==1
            if isnumeric(MarkerTable.Phenotype)
                ColName(nIP2)=strcat("IP",num2str(MarkerClusters(nIP2)));
            else
                ColName(nIP2)=MarkerClusters(nIP2);
            end

        end


        if isnumeric(MarkerTable.Phenotype)
            MarkerClusterGenes=MarkerTable.GeneSymbol(MarkerTable.Phenotype==MarkerClusters(nIP2),:);
        else
            MarkerClusterGenes=MarkerTable.GeneSymbol(MarkerTable.Phenotype==string(MarkerClusters{nIP2}));
        end

        
        
        IpANDmarker=intersect(GcaClusterGenes,MarkerClusterGenes); % IP genes that overlap with marker genes
        
        if numel(IpANDmarker)>0
            IpNOTmarker=setdiff(GcaClusterGenes,IpANDmarker); %IP genes that do not overalp with marker genes
            MarkerNOTIp=setdiff(MarkerClusterGenes,IpANDmarker); %
    
            GeneBgr=nAllGenes-numel(MarkerClusterGenes)-numel(IpNOTmarker);
    
            ContTable=table([numel(IpANDmarker);numel(IpNOTmarker)],[numel(MarkerNOTIp);GeneBgr],'VariableNames',{'InGcaCluster','NotInGcaCluster'},'RowNames',{'InMarkerCluster','NotInMarkerCluster'});
    
            
            [h,p]=fishertest(ContTable);
            
            Hmat(nIP,nIP2)=h;
            Pmat(nIP,nIP2)=p;
            Hits(nIP,nIP2)=join(IpANDmarker,',');

        else
            Hmat(nIP,nIP2)=0;
            Pmat(nIP,nIP2)=NaN;
            Hits(nIP,nIP2)={[]};
        end

    end
end

Htable=array2table(Hmat,"RowNames",RowName,"VariableNames",ColName);
Ptable=array2table(Pmat,"RowNames",RowName,"VariableNames",ColName);
Hitstable=array2table(Hits,"RowNames",RowName,"VariableNames",ColName);

figure,
hmo=heatmap(log10(Pmat),'Colormap',flipud(Default));
% title=addTitle(hmo,name,'FontSize',12);
% addXLabel(hmo,'Markers','FontSize',12);
% addYLabel(hmo,'sGCA Clusters','FontSize',12);

ax=gca;
ax.FontSize=10;
ax.ColorScaling='scaled';
ax.YDisplayLabels=RowName;
ax.XDisplayLabels=ColName;
grid off

[Prow Pcol] = find(Ptable{:,:}<0.05);

variable_names_types = [["Association", "string"];["Pval", "double"];["EnrichedGenes", "string"]];
InterTable = table('Size',[numel(Prow),size(variable_names_types,1)],'VariableNames', variable_names_types(:,1),'VariableTypes', variable_names_types(:,2));

for i=1:numel(Prow);
    InterTable.Association(i)=strjoin([ Ptable.Properties.RowNames(Prow(i)),  Ptable.Properties.VariableNames(Pcol(i))]);
    InterTable.Pval(i)=Pmat(Prow(i),Pcol(i));
    InterTable.EnrichedGenes(i)=table2array(Hitstable(Prow(i),Pcol(i)));
end

InterTable=sortrows(InterTable,"Pval",'ascend');

%save data

Info.SessionName=name;
Info.MarkerSet=fileMarker;
Info.GCAset=fileGCA;
Info.BackgroundGenes=nAllGenes;
InfoTable=struct2table(Info);

writetable(InterTable,[name,'.xlsx'],'sheet','SetAssociations');
writetable(InfoTable,[name,'.xlsx'],'sheet','Info','WriteMode','append');
saveas(ax,[name,'_sGEAmap.pdf']) %save sGEA plot as pdf
save(name);

end