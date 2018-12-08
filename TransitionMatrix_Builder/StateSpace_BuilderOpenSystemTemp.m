function [S] = StateSpace_Builder(Model, BoundaryCondition, Minimal_StateSpace)
%% State space Builder
%     This function generates the state space for the chemical master
%     equation using stoichiometry of the system

%% Check for Open species
Moiety_Conservation = sbioconsmoiety(Model,'semipos');
Open_Species_Index = find(sum(Moiety_Conservation,1) == 0);
if ~isempty(Open_Species_Index)
    if isempty(BoundaryCondition) || any(BoundaryCondition(Open_Species_Index)) == 0        
        warning('Open species exists. Define upper bound for open-species.');
        warning('Usage: Chnage the line of boundary condition in the model file with the following line');
        fprintf(strcat('Boundary_condition = zeros(1,length(CMEModel.Species));\n\n'));
        fprintf(strcat('Boundary_condition([',num2str(Open_Species_Index'),']) = ["Add a nonzero upper bound of open species as vector here"];\n\n'));
        error('Upper bound of open species is not defined');
    end
end

%% Removing non-active species
WARNING('OFF');
for i = 1:length(Model.Species)    
    delete(CMEModel.Species(i))
    BoundaryCondition(i) = [];
    fprintf('non-active species exists');
end
WARNING('ON');

%%
number_species = length(Model.Species);
number_reactions = length(Model.Reactions);
Species_Index = 1:number_species;
Reaction_Index = 1:number_reactions;

%% Setting Moiety Conservation Law.
% Converting model to a closed network by replacing null reactants/products
for i = Reaction_Index
    if isempty(get(Model.Reactions(i),'Reactants'))
         addreactant(Model.Reactions(i),['EmptyReactant_',num2str(i)]);
         OpenSpeciesFlag = 1;
    end
    if isempty(get(Model.Reactions(i),'Products'))
         addreactant(Model.Reactions(i),['EmptyProduct_',num2str(i)]);
         OpenSpeciesFlag = 1;
    end
end
Moiety_Conservation = sbioconsmoiety(Model,'semipos');
% Moiety_Conservation = Moiety_Conservation(Species_Index,:);

if ~isempty(find(Moiety_Conservation < 0, 1))
    error('Moeity conservation has negative elements')
end
if isempty(Moiety_Conservation)
    error('Moeity conservation is empty')
end

%% Initial molecular population VECTOR for generating a minimal CME
Initial_Molecular_Population = 0*Species_Index;
if Minimal_StateSpace == 1
    % In a biochemical system depending upon the initial molecular population,
    % dimension of the system varies. However, we can generate a mimal CME so
    % that each species will take part in at least one of the reaction
    % using a smaller amount of initial molecular population.
    % The floowing lines does this job by idetifying the species which are
    % neccessary for this to happen.
    
    [~, Neccessary_Species] = rref(Moiety_Conservation);
    % fprintf('To Produce a minimal CME use small values for species %d\n',Neccessary_Species);
    Initial_Molecular_Population(Neccessary_Species) = max(Moiety_Conservation,[],2);
else    
    for i = Species_Index
        Initial_Molecular_Population(i) = Model.Species(i).InitialAmount;
    end
end

% %% Non-active species
% Moiety_Conservation_Null = sbioconsmoiety(Model,'semipos');
% [RowId, ColId] = find(Moiety_Conservation_Null > 0);
% NullSpecies_Index = RowId(intersect(find(histc(RowId,unique(RowId)) == 1),find(histc(ColId,unique(ColId)) == 1)));
% 
% %% Active species
% % Open_Species_Index = find(sum(Moiety_Conservation,1) == 0);
% % Closed_Species_Index = setdiff(Species_Index',[NullSpecies_Index; Open_Species_Index]);
% Closed_Species_Index = setdiff(Species_Index',NullSpecies_Index);

%% All possible combinations of species
%     The following function identifies al possible states the system can reach
%     by means of Moiety_Conservation and upper bound.

% Closed_Moiety_Conservation = Moiety_Conservation(:,Closed_Species_Index);
% Closed_Initial_Molecular_Population = Initial_Molecular_Population(Closed_Species_Index);
% S = Combination_Finder(Closed_Moiety_Conservation, Closed_Initial_Molecular_Population);

% %% Concatinating Open species
% if ~isempty(Open_Species_Index)
%     if isempty(BoundaryCondition) || any(BoundaryCondition(Open_Species_Index)) == 0
%         error('Open Species exists')
%         warning('Open species exists. Define upper bound for open-species.');
%         warning('Usage: Chnage the commented line of boundary condition in the model file with the following line');
%         fprintf(strcat('Boundary_condition([',num2str(Open_Species_Index'),']) = ["Add a nonzero upper bound of open species as vector here"];\n\n'));
%         error('Upper bound of open species is not defined');
%     end
%     Open_S = 0:BoundaryCondition(Open_Species_Index(1));
%     for i=2:length(Open_Species_Index)
%         Open_S = combvec(Open_S,0:BoundaryCondition(Open_Species_Index(i)));
%     end
%     S = combvec(S,Open_S);
% end

% %% Reindexing back to original order as species arranged
% nstates = size(S,2);
% S(:,[Closed_Species_Index; Open_Species_Index ; NullSpecies_Index]) ...
%     = [S', repmat(Initial_Molecular_Population(NullSpecies_Index),nstates,1)];

%% Upper bound of molecular population of each species
% for a closed system, depending upon the conservation law and 
% initial molecular population, there is an upper bound on each species. 
% In open systems, we need to define the upper bound explicitely for
% generating a finite CME. For closed system algorithm  will automatically
% define the upper bopund of each species.

Conservation_Sum = Moiety_Conservation*Initial_Molecular_Population';
UB = floor(min(repmat(Conservation_Sum,1,size(Moiety_Conservation,2))./Moiety_Conservation,[],1));

%% Re-indexing: 
% While using combvec, we have to keep the vectors generated using combvec
% less exploding. In every iteration, we will generate a set 
% of column vectors using combvec and then remove some rows which are 
% not satisfying the Moeity_Conservation laws. If the first vectors are large 
% enough, then by iteratively removing the rows using Moeity_Conservation, 
% the newer vectors formed with old combvec output will be less expoloding. 
% The following re-indexing will help it.

[~,m] = sort(UB,'descend');
Moeity_Conservation_temp = Moiety_Conservation(:,m);
[~, Dependent_Species_Index] = rref(Moeity_Conservation_temp);
Independent_Species_Index = setdiff(1:size(Moeity_Conservation_temp,2),Dependent_Species_Index);

Dependent_Species_Index = m(Dependent_Species_Index);
Independent_Species_Index = m(Independent_Species_Index);

%% finding partial set of State-space
% [number_Closed_Species, number_Conseravartions] = size(Moeity_Conservation);
Independent_S = 0:UB(Independent_Species_Index(1));
for i = 2:length(Independent_Species_Index)
    Independent_S = combvec(Independent_S,0:UB(Independent_Species_Index(i)));
    
    %Removal of states that not satisfying conservation law (older/new versions of matlab)
    K = bsxfun(@minus,Conservation_Sum,Moiety_Conservation(:,Independent_Species_Index(1:i))*Independent_S);
    %   K = Conservation_Sum - A'*Moeity_Conservation(dim_Cons_Cols+1:i,:);
    Independent_S(:,any(K<0,1)) = [];    
end

%% Finding linearly dependent k species combination by solving linear problem
% Dependent_S = linsolve(Moiety_Conservation(:,Dependent_Species_Index),...
%                        bsxfun(@minus, Conservation_Sum,...
%                               Moiety_Conservation(:,Independent_Species_Index)*Independent_S));
K(:,any(K<0,1)) = [];
Dependent_S = linsolve(Moiety_Conservation(:,Dependent_Species_Index),K);
S([Dependent_Species_Index, Independent_Species_Index],:) = [Dependent_S; Independent_S];

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