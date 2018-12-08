function [D, State_Transition_Index_Matrix, SSA_propensity_matrix] = TransitionMatrix_Builder(Model, S)
% tic; 

%% Possible combination of interactions between species in all reactions
%(parallel pooling is implemented and can be used when number of reeaction is large)
[Combination_OffDiag, Combination_Diag] = Combination_reaction_Finder(Model, S);

%% Identifying the possible transitions from each state due to all reactions
State_Transition_Index_Matrix = State_Transition_Matrix_Finder(Model, S);

%% Generating off-diagonal propensitites and index
% Vecorize Reaction rates
nstates = size(S,1);
number_reactions = length(Model.Reactions);
ReactionRates = zeros(nstates,number_reactions);
for i = 1:number_reactions
    ReactionRates(:,i) = str2double(Model.Reactions(i).ReactionRate);
end

% ReactionRates(:,1) = 0.2 + 4./(1+S(:,2).^3);
% ReactionRates(:,2) = 1.09;
% ReactionRates(:,3) = 0.2 + 4./(1+S(:,1).^3);
% ReactionRates(:,4) = 1;

% Preallocation
Reaction_propensity = zeros(nstates*(1+number_reactions),1);
Reaction_propensity_index = zeros(nstates*(1+number_reactions),2);
SSA_propensity_matrix = zeros(nstates,number_reactions);
j = 0;
for i = 1:number_reactions
    % Index of possible states from where ith reaction can occur
    State_index = find(State_Transition_Index_Matrix(:,i) ~= 0);
    % Index of states after ith reaction occured
    Transition_state_index = State_Transition_Index_Matrix(State_index,i);

    l = length(Transition_state_index);
    Reaction_propensity_index(j+1:j+l,:) = [Transition_state_index, State_index];    
    Reaction_propensity(j+1:j+l) = bsxfun(@times,ReactionRates(Transition_state_index,i),...
                                                 Combination_OffDiag(State_index,i));
    
    SSA_propensity_matrix(State_index,i) = Reaction_propensity(j+1:j+l);
    
    j = j + l;
end

%% Generating diagonal propensitites and index to off diagonal propensities
Reaction_propensity_index(j+1:j+nstates,:) = [1:nstates; 1:nstates]';
Reaction_propensity(j+1:j+nstates) = -dot(Combination_Diag', ReactionRates');
j = j + nstates;

%% Removing extra preallocated area
Reaction_propensity_index(j+1:end,:) = [];
Reaction_propensity(j+1:end) = [];

%% Generating Dynamics matrix
D = sparse(Reaction_propensity_index(:,1),Reaction_propensity_index(:,2),...
           Reaction_propensity,nstates,nstates);
% nonzeros_in_D = length(Reaction_propensity);
% fprintf('Dynamics matrix generated in %.2f seconds\n\n',toc)

end