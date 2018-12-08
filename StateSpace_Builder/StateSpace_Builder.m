function [S] = StateSpace_Builder(Model, BoundaryCondition, Minimal_StateSpace)
%% State space Builder
%     This function generates the state space for the chemical master
%     equation using stoichiometry of the system

%% Check for Open species
Moiety_Conservation = sbioconsmoiety(Model,'semipos');
Open_Species_Index = find(sum(Moiety_Conservation,1) == 0);
if ~isempty(Open_Species_Index)
    if ~exist('BoundaryCondition','var') || isempty(BoundaryCondition) || any(BoundaryCondition(Open_Species_Index)) == 0        
        warning('Open species exists. Define upper bound for open-species.');
        warning('Usage: Chnage the line of boundary condition in the model file with the following line');
        fprintf(strcat('Boundary_condition = zeros(1,length(CMEModel.Species));\n\n'));
        fprintf(strcat('Boundary_condition([',num2str(Open_Species_Index'),']) = ["Add a nonzero upper bound of open species as vector here"];\n\n'));
        error('Upper bound of open species is not defined');
    end
end

%% Removing non-active species
warning('OFF', 'SimBiology:DELETE_SPECIES_BEING_USED');
for i = 1:length(Model.Species)    
    delete(Model.Species(i))
end

number_species = length(Model.Species);
%% Setting conservation Law.
%     For a closed system there are a number of conservartion laws exists. For
%     open system, there may be some but not neccessary. We identify the
%     possible conservation laws in both systems using the following function
Moiety_Conservation = sbioconsmoiety(Model,'semipos');
if ~isempty(find(Moiety_Conservation < 0, 1))
    error('Moeity conservation has negative elements')
end

%% Initial molecular population for generating a minimal CME
%     In a biochemical system depending upon the initial molecular population,
%     dimension of the system varies. However, we can generate a mimal CME so
%     that each species will take part in at least one of the reaction
%     using a smaller amount of initial molecular population.
%     The floowing lines does this job by idetifying the species which are
%     neccessary for this to happen.

[~, Neccessary_Species] = rref(Moiety_Conservation);
% fprintf('To Produce a minimal CME use small values for species %d\n',Neccessary_Species);
Minimal_Initial_Molecular_Population = zeros(1,number_species);
Minimal_Initial_Molecular_Population(Neccessary_Species) = max(Moiety_Conservation,[],2);

%% Vectorize Initial molecular population
if Minimal_StateSpace == 1
    Initial_Molecular_Population = Minimal_Initial_Molecular_Population;
else
    Initial_Molecular_Population = zeros(1,number_species);
    for i = 1:number_species
        Initial_Molecular_Population(i) = Model.Species(i).InitialAmount;
    end
end

%% Removing non-active species and correspinding relations
% [RowId, ColId] = find(Moiety_Conservation > 0);
% NullSpecies_Index = RowId(intersect(find(histc(RowId,unique(RowId)) == 1),find(histc(ColId,unique(ColId)) == 1)));
Open_Species_Index = find(sum(Moiety_Conservation,1) == 0);
Closed_Species_Index = setdiff((1:number_species)', Open_Species_Index);

%% All possible combinations of species
%     The following function identifies al possible states the system can reach
%     by means of Moiety_Conservation and upper bound.

Closed_Moiety_Conservation = Moiety_Conservation(:,Closed_Species_Index);
Closed_Initial_Molecular_Population = Initial_Molecular_Population(Closed_Species_Index);
Conservation_Sum = Moiety_Conservation*Initial_Molecular_Population';
S = IntegerSolutionsAll(Closed_Moiety_Conservation, Conservation_Sum);

%% Concatinating Open species
if ~isempty(Open_Species_Index)
    Open_S = 0:BoundaryCondition(Open_Species_Index(1));
    for i=2:length(Open_Species_Index)
        Open_S = combvec(Open_S,0:BoundaryCondition(Open_Species_Index(i)));
    end
    S = combvec(S,Open_S);
end

S = S';
%% Reindexing back to original order as species arranged
S(:,[Closed_Species_Index; Open_Species_Index]) = S;

%% Check for absence of species
%     checking for zero column in S so that each species is present
%     in atleast one state

Zero_species = find(sum(S) == 0);
if ~isempty(Zero_species)
    display(Moiety_Conservation);
    display(Initial_Molecular_Population);
    warning('Not all species states found.')
    warning('Check negative values in Conservation matrix');
    warning('Check whether intial population has all necessary species');
    error('Species %d are not present in state space\n',Zero_species)
end
end