clear all
close all
clc
%% add path to toolbox and intialize
addpath('/Users/bpb/matlabtoolbox/plotting_functions/plot2svg_20101128/plot2svg_20101128/')
addpath('/Users/bpb/matlabtoolbox/metabolomics_IO_functions/xml_io_tools_2010_05_04/')
addpath('/Users/bpb/matlabtoolbox/math_models/cobra/')
addpath('/Users/bpb/matlabtoolbox/metabolomics_IO_functions/readtext/')
initCobraToolbox
%% open the sbml file and replace the & with "and" then parse using xml_read
fid = fopen('ModifiedShewyMR1_Mar2009 (3).xml','r');
str=char(fread(fid,'*uchar')');
fclose(fid)
str=regexprep(str,'&','and');
fid=fopen('temp.xml','w');
fprintf(fid,'%s',str);
fclose(fid);
str=xml_read('temp.xml');
%%
newmodel=MyParseSBMLtoCobra(str);
%% make rxnGeneMap
% get unique genes that don't participate in OR reactions.  We'll deal with
% these later.
tic
c=1;
genes=cell(1e6,1);
for i = 1:length(newmodel.grRules);
    if isempty(strfind(newmodel.grRules{i},' or '))
        a=regexp(newmodel.grRules{i},'SO\d+','match');
        for j = 1:length(a)
            genes{c}=a{j};
            c=c+1;
        end
    end
end
genes(c:end)=[];
newmodel.genes=unique(genes);
newmodel.rxnGeneMat=zeros(size(newmodel.rxns,1),length(newmodel.genes));
for i = 1:length(newmodel.grRules);
    if isempty(strfind(newmodel.grRules{i},' or '))
        a=regexp(newmodel.grRules{i},'SO\d+','match');
        [u ua ub]=intersect(a,newmodel.genes);
        newmodel.rxnGeneMat(i,ub)=1;
    end
end
toc
%% remove all the bulk metabolites
r=regexp(newmodel.mets,'_b$','match');
xx=zeros(size(r));
for i = 1:length(xx)
    if ~isempty(r{i})
        xx(i)=1;
    end
end
% xx is the indices of metabolites that need to be removed
newmodel.mets(xx==1)=[];
newmodel.metNames(xx==1)=[];
newmodel.metCharge(xx==1)=[];
newmodel.S(xx==1,:)=[];

%%  Add growth dependant ATP maintenance requirements
xx=strmatch('SO_BIOMASSMACRO_DM_NOATP2',newmodel.rxns);
newmodel.c(xx)=1;
GAM=220.22;
spe={'ATP_C','ADP_C','H2O_C','PI_C','H_C'};
ua=zeros(size(spe));
for i = 1:length(spe)
    ua(i)=strmatch(spe{i},upper(newmodel.mets));
end
edges=[-GAM GAM -GAM GAM GAM];
newmodel.S(ua,xx)=edges;
%% Change some formatting so that it is like the other published Shewy Model
newmodel.mets=regexprep(newmodel.mets,'_c$','[c]');
newmodel.mets=regexprep(newmodel.mets,'_e$','[e]');
newmodel.mets=upper(regexprep(newmodel.mets,{'formate','4_2','p-dd'},{'for','4:2','p-D,D'}))

%% make a quick printout of the reactions
clc
xx=strmatch('EX_',newmodel.rxns);
fid=fopen('temp.tab','w');
fprintf(fid,'Reaction Names\tReaction ID\n');
for i = 1:length(xx)
    fprintf(fid,'%s\t%s\n',newmodel.rxnNames{xx(i)},newmodel.rxns{xx(i)});
end
fclose('all')
%% save model as a backup to a matlab structure
save('ShewyModel','newmodel');
%% setup the loader of all the media compositions
r=readtext('MR1_experiments_withgrowthrate.tab','\t');
%%
Conditions.Chip=[r{2:end,1}];
Conditions.InUse=[r{2:end,5}];
Conditions.Group=r(2:end,2);
Conditions.GrowthMethod=r(2:end,3);
Conditions.Info=r(2:end,4);
Conditions.CarbonSource=r(2:end,6);
Conditions.AllowedFlux=r(2:end,7:end);
for i = 1:length(Conditions.Info)
    Conditions.Info{i}=sprintf('%s: %s',Conditions.Group{i},Conditions.Info{i});
end
%% read in the gene data
fid = fopen('MR1_genes.tab','r');
str=fread(fid);
fclose(fid);
str=char(str)';
size(str);
a=textscan(str,'%s','delimiter','\n');
fn=textscan(a{1}{1},'%s','delimiter','\t');
fn=fn{1};
for i= 1:length(a{1})-1
    temp=textscan(a{1}{i+1},'%s','delimiter','\t');
    temp=temp{1};
    for j = 1:length(temp)
        geneinfo.(fn{j}){i}=temp{j};
    end
end
%% spin out all the models
xx=strmatch('EX_',newmodel.rxns);
newmodel.lb(xx)=0;
newmodel.ub(xx)=1000; % zero out all exchange reactions
spe={'EX_CA2_E',...
    'EX_CL_E',...
    'EX_H2O_E',...
    'EX_MG2_E',...
    'EX_H_E',...
    'EX_CO2_E',...
    'EX_NA1_E',...
    'EX_K_E'}

for i = 1:length(spe)
    xx=strmatch(spe{i},newmodel.rxns);
    newmodel.lb(xx)=-1000;
end

xx5=find(Conditions.InUse~=0);
g=zeros(size(Conditions.InUse));
for i = 1:length(xx5)
    tempmodel=newmodel;
    xx=strmatch(Conditions.CarbonSource{xx5(i)},newmodel.rxns);
    tempmodel.lb(xx)=-10;
    for j = 1:size(Conditions.AllowedFlux,2)
        if ~isempty(Conditions.AllowedFlux{xx5(i),j})
            xx=strmatch(Conditions.AllowedFlux{xx5(i),j},newmodel.rxns);
            tempmodel.lb(xx)=-1000;
        end
    end
    temp=optimizeCbModel(tempmodel);
    g(xx5(i))=temp.f;
end
%% compare to chris's D,L lactate - stationary (high Calcium)  growth
idx=strncmp('carbon source: D,L lactate - stationary (high Calcium) ',Conditions.Info,25)
g(idx)
%% Run model-wide single gene deletion analysis using FBA:
% make sure that when a gene goes with multiple reactions, both reactions
% are cut for gene deletion simulation.  Keep the reaction cutter
% gene-oindependent.s
xx=strmatch('EX_',newmodel.rxns);
newmodel.lb(xx)=0;
newmodel.ub(xx)=1000; % zero out all exchange reactions
spe={'EX_CA2_E',...
    'EX_CL_E',...
    'EX_H2O_E',...
    'EX_MG2_E',...
    'EX_H_E',...
    'EX_CO2_E',...
    'EX_NA1_E',...
    'EX_K_E'}
% waste ammonium
for i = 1:length(spe)
    xx=strmatch(spe{i},newmodel.rxns);
    newmodel.lb(xx)=-1000;
    newmodel.ub(xx)=1000;
end

xx5=find(Conditions.InUse~=0);
myflux=zeros(length(xx5),length(tempmodel.genes));
myrflux=zeros(length(xx5),length(tempmodel.rxns));

for ii = 1:length(xx5)
    tempmodel=newmodel;
    xx=strmatch(Conditions.CarbonSource{xx5(ii)},newmodel.rxns);
    tempmodel.lb(xx)=-10;
    for j = 1:size(Conditions.AllowedFlux,2)
        if ~isempty(Conditions.AllowedFlux{xx5(ii),j})
            xx=strmatch(Conditions.AllowedFlux{xx5(ii),j},newmodel.rxns);
            tempmodel.lb(xx)=-1000;
        end
    end
    for iii = 1:length(tempmodel.genes)
        xx=find(tempmodel.rxnGeneMat(:,iii)~=0);
        tempmodel2=tempmodel;
        tempmodel2.lb(xx)=-1e-9;
        tempmodel2.ub(xx)=1e-9;
        s=optimizeCbModel(tempmodel2);
        myflux(ii,iii)=s.f;
%         myrflux(ii,xx)=s.f;
        [iii ii]
    end
end

for ii = 1:length(xx5)
    tempmodel=newmodel;
    xx=strmatch(Conditions.CarbonSource{xx5(ii)},newmodel.rxns);
    tempmodel.lb(xx)=-10;
    for j = 1:size(Conditions.AllowedFlux,2)
        if ~isempty(Conditions.AllowedFlux{xx5(ii),j})
            xx=strmatch(Conditions.AllowedFlux{xx5(ii),j},newmodel.rxns);
            tempmodel.lb(xx)=-1000;
        end
    end
    for iii = 1:length(newmodel.rxns)
        tempmodel2=tempmodel;
        tempmodel2.lb(iii)=-0.0001;
        tempmodel2.ub(iii)=0.0001;
        s=optimizeCbModel(tempmodel2);
        myrflux(ii,iii)=s.f;
        [iii ii]
    end
end
save KOfluxMediaMatrix myflux myrflux xx5


%% export table of knockouts and knockreactionsout as a function of growth conditions
fid = fopen('KO_predicted_growth2.tab','w');
fprintf(fid,'ChipID\t');
fprintf(fid,'%s\t',newmodel.genes{:});
fprintf(fid,'\n');
for i = 1:length(xx5)
    fprintf(fid,'%d\t',Conditions.Chip(xx5(i)));
    fprintf(fid,'%d\t',myflux(i,:));
    fprintf(fid,'\n');
end
fclose('all');

fid = fopen('KO_reactions_predicted_growth2.tab','w');
fprintf(fid,'ChipID\t');
fprintf(fid,'%s\t',newmodel.rxnNames{:});
fprintf(fid,'\n\t');
fprintf(fid,'%s\t',newmodel.rxns{:});
fprintf(fid,'\n');
for i = 1:length(xx5)
    fprintf(fid,'%d\t',Conditions.Chip(xx5(i)));
    fprintf(fid,'%d\t',myrflux(i,:));
    fprintf(fid,'\n');
end
fclose('all');
% make a quick printout of the reactions
fid=fopen('ModelReactions2.tab','w');
fprintf(fid,'Reaction Names\tReaction ID\tReaction\tGeneMapStr\tGenes\n');
for i = 1:size(newmodel.S,2)
    fprintf(fid,'%s\t%s\t',newmodel.rxnNames{i},newmodel.rxns{i});
    xx=find(newmodel.S(:,i)<0);
    for j = 1:length(xx)  
        if abs(newmodel.S(xx(j),i))~=1
        fprintf(fid,'%d %s',abs(newmodel.S(xx(j),i)),newmodel.mets{xx(j)});
        else
        fprintf(fid,'%s',newmodel.mets{xx(j)});
        end            
        if j ~=length(xx)
            fprintf(fid,' + ');
        else
            fprintf(fid,'->');
        end
    end
    xx=find(newmodel.S(:,i)>0);
    for j = 1:length(xx)  
if abs(newmodel.S(xx(j),i))~=1
        fprintf(fid,'%d %s',abs(newmodel.S(xx(j),i)),newmodel.mets{xx(j)});
        else
        fprintf(fid,'%s',newmodel.mets{xx(j)});
end
if j ~=length(xx)
            fprintf(fid,' + ');
        else
            fprintf(fid,'->');
        end
    end
    fprintf(fid,'\t');
    fprintf(fid,'%s\t',newmodel.grRules{i});
    xx=find(newmodel.rxnGeneMat(i,:)~=0);
    fprintf(fid,'%s;',newmodel.genes{xx});
    fprintf(fid,'\n');
end
fclose('all')
%%
try
    load mr2fitness
catch
r=readtext('MR1_282_data.tab','\t');
fitness.gi=[r{2:end,1}]
fitness.genes=r(2:end,3)
fitness.experiment=[r{1,4:end}]
fitness.value=cell2mat(r(2:end,4:end));

r=readtext('MR1_282_metadata.tab','\t');
 [u ua ub]=intersect(fitness.experiment,[r{2:end,1}]);
    fitness.experiment=fitness.experiment(ua);
    fitness.group=r(2:end,2);
    fitness.group=fitness.group(ub);
    fitness.info=r(2:end,3);
    fitness.info=fitness.info(ub);
    save mr2fitness fitness
end

%% map the genes to reaction to calculate a reaction specific delta-fitness
load mr2fitness
load ShewyModel
myrfitness=zeros(length(newmodel.rxns),length(fitness.experiment));
for i = 1:length(newmodel.rxns)
    xx=find(newmodel.rxnGeneMat(i,:)~=0);
    temp=zeros(length(xx),length(fitness.experiment));
    [u ua ub]=intersect(newmodel.genes(xx),fitness.genes);   
        myrfitness(i,:)=mean(fitness.value(ub,:),1);
        i
end
%% play with flux matrix
load KOfluxMediaMatrix
load ShewyModel
size(myrflux)
s1=find(sum(myrflux,1)==0)
s1=find(sum(myrflux,2)~=0)
costr=Conditions.Info(xx5(s1));
chipused=Conditions.Chip(xx5(s1));
temp=myrflux(s1,:);
tempfit=myrfitness(:,s1)';
% temp is the flux matrix [94 media x 573 genes]


% u=unique(newmodel.subSystems) % plot clustergram of simulations for tca cycle
% xx=strmatch('SUBSYSTEM: Glycolysis/Gluconeogenesis',newmodel.subSystems);
% xx=cat(1,xx,strmatch('SUBSYSTEM: Citrate Cycle (TCA)',newmodel.subSystems));
xx1=strmatch('SUBSYSTEM: Citrate Cycle (TCA)',newmodel.subSystems);
xx2=strmatch('SUBSYSTEM: Glycolysis/Gluconeogenesis',newmodel.subSystems);
xx1=union(xx1,xx2);
xx2=strmatch('SUBSYSTEM: Glycine and Serine Metabolism',newmodel.subSystems);
xx1=union(xx1,xx2);
xx2=strmatch('SUBSYSTEM: Valine, Leucine, and Isoleucine Metabolism',newmodel.subSystems);
xx=union(xx1,xx2);
% xx contains the list of reactions that are for a
r=temp(:,xx);
rf=tempfit(:,xx);
%m=find(sum(newmodel.rxnGeneMat(xx,:),1)~=0);
%r=temp(:,m)
% get unique medias
[u ua] = unique(r,'rows');
g=newmodel.rxnNames(xx);%  genes(m);
cg=clustergram(r(ua,:)','rowlabels',g,'columnlabels',costr(ua))
z=r(ua,:)';
zfit=rf(ua,:)';
x1=get(cg,'ColumnLabels')
y1=get(cg,'RowLabels')


%% try and make a pretty surface plot with better labels.
temp2=r(ua,:)';
temp2(temp2>1)=1;
ux1=zeros(size(x1));
uy1=zeros(size(y1));
for i = 1:length(x1);
    t=strmatch(x1{i},costr(ua),'exact');
    ux1(i)=t(1);
end
for i = 1:length(y1)
    uy1(i)=strmatch(y1{i},g,'exact')
end
temp=temp2(uy1,ux1);
tempfit=zfit(uy1,ux1);
tempfit(~isfinite(tempfit))=0;%min(tempfit(:));
% temp(temp>2)=2;
% temp=temp-min(temp(:));
% temp=temp/max(temp(:));
MakeNiceClustergramPlot(temp,x1,y1)
% caxis([0 2.5])
c=colormap(jet);
c(1,:)=[0 0 0];
colormap(c)
figure
MakeNiceClustergramPlot(tempfit,x1,y1)
c=colormap(jet);
c(1,:)=[0 0 0];
colormap(c)
figure
% tf2=temp;
% tf2=tf2>0.1;
tf2=tempfit;
tf2(abs(temp)<0.1)=min(tf2(:));
% tf2(tf2<-1)=-1;
% tf2=tf2+1;
MakeNiceClustergramPlot(tf2,x1,y1)
c=colormap(jet);
c(1,:)=[1 1 1];
colormap(c)


% plot the fi
% use plot2svg to save it
xx=find(temp(:)>0.01 & tempfit(:)~=0);
figure, subplot(2,1,1), hist(tempfit(xx),100)
xx=find(temp(:)<0.01 & tempfit(:)~=0);
subplot(2,1,2), hist(tempfit(xx),100)



%% make a small graph of the flux for a particular pathway
xx1=strmatch('SUBSYSTEM: Citrate Cycle (TCA)',newmodel.subSystems);
length(xx1)
xx2=strmatch('SUBSYSTEM: Glycolysis/Gluconeogenesis',newmodel.subSystems);

xx1=union(xx1,xx2);
length(xx1)
xx2=strmatch('SUBSYSTEM: Glycine and Serine Metabolism',newmodel.subSystems);
xx1=union(xx1,xx2);
xx2=strmatch('SUBSYSTEM: Valine, Leucine, and Isoleucine Metabolism',newmodel.subSystems);
xx1=union(xx1,xx2);


length(xx1)
% xx1=1:length(newmodel.rxns);
rstr=newmodel.rxns(xx1);
rstr2=newmodel.rxnNames(xx1);
X=newmodel.S(:,xx1);
flux=randn(size(xx1));%solution.x(xx1);
xx2=find(max(abs(X),[],2)~=0);
X=X(xx2,:);
mstr=newmodel.mets(xx2);
mstr2=newmodel.metNames(xx2);
% make copies of C02, H2O, NAD, NADH, etc
cstr={'NH4[C]','CO2[C]','ACCOA[C]','COA[C]','NAD[C]','NADH[C]','H[C]','H2O[C]','PI[C]','ATP[C]','ADP[C]'};
[u ua ub]=intersect(cstr,mstr);
u2=mstr2(ub);
temp=X(ub,:);
X(ub,:)=[];
mstr(ub)=[];
mstr2(ub)=[];
for i = 1:size(temp,1)
    xx=find(temp(i,:)~=0);
    for j = 1:length(xx)
        t1=zeros(1,size(temp,2));
        t1(xx(j))=temp(i,xx(j));
        X=cat(1,X,t1);
        mstr=cat(1,mstr,[u{i}]);
        mstr2=cat(1,mstr2,[u2{i}]);
    end
end

% turn the stoichiometric matrix into a connectivity matrix
str=cat(1,rstr,mstr);
namestr=cat(1,rstr2,mstr2);
cm=zeros(length(str));
% find each metabolite that feeds into a reaction (ie: neg stoich)
[ii jj]=find(X<0);
ii=ii+length(rstr);
for i = 1:length(ii)
    cm(ii(i),jj(i))=1;
end
% find each metabolite that feeds out of a reaction (ie: pos stoich)
[ii jj]=find(X>0);
ii=ii+length(rstr);
for i = 1:length(ii)
    cm(jj(i),ii(i))=1;
end

bg=biograph(cm)
for i = 1:length(bg.Nodes)
    bg.Nodes(i).Label=str{i};
    bg.Nodes(i).Description=namestr{i};
    if i>length(rstr)
        bg.Nodes(i).Shape='box';
        bg.Nodes(i).Color=[0.1 1 0.1];
        bg.Nodes(i).LineColor=[0 0 0];
    else
        bg.Nodes(i).Shape='box';
        bg.Nodes(i).Color=[1 0.1 0.1];
        bg.Nodes(i).LineColor=[0 0 0];
        bg.Nodes(i).FontSize=12;
    end
end
% set the edges color relative to the flux
figure(1)
map=colormap(jet); % goes red to yellow
close all
idx1=linspace(0,max(abs(flux)),size(map,1));
% idx1=idx1(end:-1:1);
% figure(1)
% map2=colormap(summer); % goes yellow to green
% close all
% idx2=linspace(min(flux),0,size(map,1));
% idx2=idx2(end:-1:1);
for i = 1:length(bg.Edges)
    temp=bg.Edges(i).ID;
    a=regexp(temp,'Node \d+','match');
    a=regexprep(a,'Node ','');
    a=str2double(a);
    xx=find(a<length(rstr));
    if ~isempty(xx)
        [m mx]=min(abs(abs(flux(a(xx)))-idx1));
        
        if flux(a(xx))>0
            bg.Edges(i).Weight=4;
            bg.Edges(i).LineColor=map1(mx,:);
        elseif flux(a(xx))<0
            bg.Edges(i).Weight=4;
            
            bg.Edges(i).LineColor=map1(mx,:);
        elseif flux(a(xx))==0
            bg.Edges(i).Weight=4;
            
            bg.Edges(i).LineColor=[0.5 0.5 0.5];
            
        end
    end
    
end


view(bg)


%%

    filename=['C:\matlabtoolbox\math_models\models\ShewanellaOdeinsis\',...
        'Shewanella_odeinsis_Pinchuck_Sup_PlosCompBio_FixedBiomass.xls']
    %     [n,t,r]=xlsread(filename,'reactions')
    model = xls2model(filename);%,biomassRxnEquation)
    xx=strmatch('Biomass',model.rxns,'exact')
    model.c(xx)=1;
    solution=optimizeCbModel(model)
    save mr2model model


%% thoughts
% don't overthink
% interested in genes where knockout is predicted to grow, but fitness is
% clearly (ie < -2) inhibit growth or vice-versa (fitess > -0.5) but the
% model predicts no growth

% do n-acetyl glucosamine and lactate as sole-carbon sources with ammonium
% as the sole nitrogen source.  Only use aerobic growth

% automating using the fitess data to improve gene-reaction annotations

% 1) Get the model to grow on all the conditions
% 2) Quantify the agreement and discrepancy between the model and the fitness
% data
% 3) Clean up how I'm knocking out genes/reactions
% 4) bootstrap effectiveness of model improvements
% 5) put the metrics on susbsytesm (ie: find subsystems where KO makes
% model not grow, but fitness shows growth on a variety of conditions).
%
% Morgan will tabulate the modifications to make to the model
%
% next, I will figure out how to run it with Morgan's changes
%
% We need to send an email with our progress to Jenny and Chris, and ask
% them if they would be willing to try the growmatch on the original and
% refined model

% debugging model code and scripts

% making the model grow in all the conditions we have

% making metrics for how well the model agrees

% editing the model to improve the fit (adding reactions, changing gene
% associations, removing reactions, adding or removing componnents to the
% biomass).  There are also regulatory and LB/UB consttrains, first simple
% stuff.

% given that fitness doesn't reflect significant growth rates on a given
% media (ie: you can have a high fitness while having a very low growth
% rate).  Only KO situations can be considered.  On a given media, this KO
% performs much better than another KO.

% take susbsytems and show how various subsystems are in good agreement or
% bad agreement.
% 1) select a subsystem and
% 2) the rows will be enzyme in that subsystem
% 3) the columns will be a "quality score" for how it matches KO data
% cluster this matrix and show a network layout for each subsystem

% for next week
% allow waste disposal
% show plots for reactions not genes.

% given that there is a unique set of predicted fluxes for the various
% medias (ie: all lactates will yield the same flux).  Are other media
% formulations identical?  For example, is there another media formulation
% that gives identical results to lactate?  How does this similarity
% compare to the similarities observed for the KO fitness data?

% the question is how to have metabolites that are present in one subsystem
% connect to another subsystem.  For example, if none connect than you will
% only have subsystem networks.  If all connect, you will eventually have a
% hairball.

