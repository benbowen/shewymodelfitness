function newmodel = MyParseSBMLtoCobra(str)
% str.model.listOfSpecies.species.ATTRIBUTE is the struct that lists the compounds and their compartment.
%
% str.model.listOfUnitDefinitions.unitDefinition.listOfUnits.unit.ATTRIBUTE
% lists the units (mole, gram, second) and appropriate multipliers
%
% str.model.listOfCompartments.compartment.ATTRIBUTE lists the two
% compartments
%
%                   rxns: {870x1 cell}
%               rxnNames: {870x1 cell}
%             subSystems: {870x1 cell}
%                     lb: [870x1 double]
%                     ub: [870x1 double]
%                    rev: [870x1 double]
%                      c: [870x1 double]
%                      b: [713x1 double]
%                      S: [713x870 double]
%             rxnGeneMat: [870x783 double]
%                  rules: {870x1 cell}
%                grRules: {870x1 cell}
%                  genes: {783x1 cell}
%           rxnECNumbers: {870x1 cell}
%               proteins: {870x1 cell}


%% make model.mets
% make it look like this
%     'uri[e]'
%     'ttdca[c]'
tic
newmodel=[]
newmodel.mets=cell(length(str.model.listOfSpecies.species(:)),1);
newmodel.metNames=cell(length(str.model.listOfSpecies.species(:)),1);
% newmodel.metFormulas=cell(length(str.model.listOfSpecies.species(:)),1);
% newmodel.metFormulasNeutral=cell(length(str.model.listOfSpecies.species(:)),1);
% newmodel.metKEGGID=cell(length(str.model.listOfSpecies.species(:)),1);
newmodel.metCharge=zeros(length(str.model.listOfSpecies.species(:)),1);

for i = 1:length(newmodel.mets)
    newmodel.mets{i}=[str.model.listOfSpecies.species(i).ATTRIBUTE.id];%,'[',str.model.listOfSpecies.species(i).ATTRIBUTE.compartment,']'];
    newmodel.metNames{i}=str.model.listOfSpecies.species(i).ATTRIBUTE.name;
    newmodel.metCharge(i)=str.model.listOfSpecies.species(i).ATTRIBUTE.charge;
end
toc
[u ua]=unique(newmodel.mets);
newmodel.mets=u;
newmodel.metNames=newmodel.metNames(ua);
newmodel.metCharge=newmodel.metCharge(ua);
% newmodel.mets=regexprep(newmodel.mets,'Extracellular','e');
% newmodel.mets=regexprep(newmodel.mets,'Cytosol','c');
% make model.rxns
%
% lists the species and stoichiometry listofproducts and listofreactants
%                      b: [713x1 double]
%                      S: [713x870 double]
%             rxnGeneMat: [870x783 double]
%                  rules: {870x1 cell}
%                grRules: {870x1 cell}
%                  genes: {783x1 cell}
%           rxnECNumbers: {870x1 cell}
%               proteins: {870x1 cell}

tic
newmodel.rxns=cell(length(str.model.listOfReactions.reaction),1);
newmodel.rxnNames=cell(length(str.model.listOfReactions.reaction),1);
newmodel.subSystems=cell(length(str.model.listOfReactions.reaction),1); % this could contain pathway classes for hte model
newmodel.rules=cell(length(str.model.listOfReactions.reaction),1);
newmodel.grRules=cell(length(str.model.listOfReactions.reaction),1);
% newmodel.rxnECNumbers=cell(length(str.model.listOfReactions.reaction),1);
newmodel.proteins=cell(length(str.model.listOfReactions.reaction),1);

newmodel.lbunits=cell(length(str.model.listOfReactions.reaction),1);
newmodel.ubunits=cell(length(str.model.listOfReactions.reaction),1);
newmodel.fluxunits=cell(length(str.model.listOfReactions.reaction),1);


newmodel.lb=zeros(length(str.model.listOfReactions.reaction),1);
newmodel.ub=zeros(length(str.model.listOfReactions.reaction),1);
newmodel.rev=zeros(length(str.model.listOfReactions.reaction),1);
newmodel.c=zeros(length(str.model.listOfReactions.reaction),1);
newmodel.flux=zeros(length(str.model.listOfReactions.reaction),1);

newmodel.S=zeros(length(newmodel.mets),length(str.model.listOfReactions.reaction));

for i = 1:length(newmodel.rxns)
    newmodel.rxns{i}=str.model.listOfReactions.reaction(i).ATTRIBUTE.id;
    newmodel.rxnNames{i}=str.model.listOfReactions.reaction(i).ATTRIBUTE.name;
    newmodel.rev(i)=strcmp(str.model.listOfReactions.reaction(i).ATTRIBUTE.reversible,'true');
    newmodel.subSystems{i}=str.model.listOfReactions.reaction(i).notes.html_COLON_p{3};
    newmodel.rules{i}=str.model.listOfReactions.reaction(i).notes.html_COLON_p{2};
    newmodel.grRules{i}=str.model.listOfReactions.reaction(i).notes.html_COLON_p{1};
    newmodel.proteins{i}=str.model.listOfReactions.reaction(i).notes.html_COLON_p{2};
    newmodel.lb(i)=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(1).ATTRIBUTE.value;
    newmodel.c(i)=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(3).ATTRIBUTE.value;
    newmodel.ub(i)=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(2).ATTRIBUTE.value;
    newmodel.flux(i)=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(4).ATTRIBUTE.value;
    
    newmodel.lbunits{i}=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(1).ATTRIBUTE.units;
    newmodel.ubunits{i}=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(2).ATTRIBUTE.units;
    newmodel.fluxunits{i}=str.model.listOfReactions.reaction(i).kineticLaw.listOfParameters.parameter(4).ATTRIBUTE.units;
    if ~isempty(str.model.listOfReactions.reaction(i).listOfReactants)
        for j = 1:length(str.model.listOfReactions.reaction(i).listOfReactants.speciesReference)
            sp=str.model.listOfReactions.reaction(i).listOfReactants.speciesReference(j).ATTRIBUTE.species;
            st=str.model.listOfReactions.reaction(i).listOfReactants.speciesReference(j).ATTRIBUTE.stoichiometry;
            xx1=strcmp(sp,newmodel.mets);
            newmodel.S(xx1,i)=st*-1;
        end
    end
    if ~isempty(str.model.listOfReactions.reaction(i).listOfProducts)
        for j = 1:length(str.model.listOfReactions.reaction(i).listOfProducts.speciesReference)
            sp=str.model.listOfReactions.reaction(i).listOfProducts.speciesReference(j).ATTRIBUTE.species;
            st=str.model.listOfReactions.reaction(i).listOfProducts.speciesReference(j).ATTRIBUTE.stoichiometry;
            xx1=strcmp(sp,newmodel.mets);
            newmodel.S(xx1,i)=st;
        end
    end
    
    
end
