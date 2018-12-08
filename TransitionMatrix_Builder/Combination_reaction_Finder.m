function [Combination_OffDiag, Combination_diag] = Combination_reaction_Finder(Model, S)
%
Stoichiometry = Reactant_stoichiometry(Model);
nstates = size(S,1);
number_reactions = length(Model.Reactions);
Combination_OffDiag = zeros(nstates,number_reactions);
Combination_diag = zeros(nstates,number_reactions);
for i = 1:number_reactions %(USE parfor for parallel pooling)
    % i-th reaction stoichiometry
    Reaction_i = Stoichiometry(:,i); % to avoid overhead in parallel pool
    % Positions of reactants
    Reactant_Position = find(Reaction_i < 0);
    % Number of molecules (Reactants) needed for reaction i
    Reactant_change = -Reaction_i(Reactant_Position)';
    % To avoids negative values in gammaln function
    S_nonneg_position = find(all(bsxfun(@minus,S(:,Reactant_Position),Reactant_change) >= 0,2));
    % Possible states that the reaction i occurs
    S_noneg_reactants = S(S_nonneg_position,Reactant_Position);
    Reactant_change = repmat(Reactant_change, length(S_nonneg_position),1);    
    % Finds nchoosek(S,R): number of ways reaction i can occur with S
    % molecule under R reactant consumption
    % gammaln function does the nchoosek function without overflow
    Combination_OffDiag(S_nonneg_position,i) = prod(round(exp(gammaln(S_noneg_reactants+1)...
        -gammaln(Reactant_change+1)-gammaln(S_noneg_reactants-Reactant_change+1))),2);
    
    % Dynamic matrix's diagonal element combinations
    Combination_diag(:,i) = prod(S(:,Reactant_Position),2);
end
end


% function [Combination_OffDiag, Combination_diag] = Combination_reaction_Finder(Model, S)
% 
% %%
% Stoichiometry = Reactant_stoichiometry(Model);
% nstates = size(S,1);
% number_reactions = length(Model.Reactions);
% Combination_OffDiag = zeros(nstates,number_reactions);
% Combination_diag = zeros(nstates,number_reactions);
% for i = 1:number_reactions %(USE parfor for parallel pooling)
%     % i-th reaction stoichiometry
%     Reaction_i = Stoichiometry(:,i); % to avoid overhead in parallel pool
%     % Positions of reactants
%     Reactant_Position = find(Reaction_i < 0);
% 
% 
%     % To generate vectors for combination products
%     Comb = zeros(nstates,length(Reactant_Position));
%     % To avoid negative values in using nchoosek
%     for j = 1:length(Reactant_Position)
%         % Number of molecules (Reactants) needed for reaction i
%         Reactant_change = -Reaction_i(Reactant_Position(j));
%         % To avoids negative values in gammaln function
%         S_non_negative_position = find(S(:,Reactant_Position(j)) - Reactant_change >= 0);
%         % Possible states that the reaction i occurs
%         S_temp = S(S_non_negative_position,Reactant_Position(j));
%         % Finds nchoosek(S,R): number of ways reaction i can occur with S
%         % molecule under R reactant consumption
%         % gammaln function does the nchoosek function without overflow
%         Comb(S_non_negative_position,j) = round(exp(gammaln(S_temp+1)...
%             -gammaln(Reactant_change+1)-gammaln(S_temp-Reactant_change+1)));
%     end
% 
%     %Possible Combinations of the occurance of reaction i from all states
%     Combination_OffDiag(:,i) = prod(Comb,2);
%     % Dynamic matrix's diagonal element combinations
%     Combination_diag(:,i) = prod(S(:,Reactant_Position),2);
% end
% end