function State_Transition_Index_Matrix = State_Transition_Matrix_Finder(S,Stoichiometry)
%% The function identifies the transitioned states from all possible states
    % due to different reactions. 
% tic
global nstates number_reactions number_species
% (States + Stoichiometry) is implememnted in a single line to avoid memory
% loss and avoid for loop whcih can also be identified through-

% S_temp = repmat(S,number_reactions,1);
% Stoichiometry_temp = kron(Stoichiometry',ones(nstates,1));
% B = S_temp + Stoichiometry_temp;

B = repmat(S,number_reactions,1) + kron(Stoichiometry',ones(nstates,1));

% Removing unrealistic states
B(any(B < 0,2),:) = 0;

% Identifying unique states to improve the ismember function
[Unique_states,~,Index_in_B] = unique(B,'rows');

% Identfying the index of unique states in S
[~,Index_in_S] = ismember(Unique_states,[S; zeros(1,number_species)],'rows');

% Idnetfying the position of each state in S after each reactons 
Index_in_S(Index_in_S == nstates+1) = 0;
Index_in_S = Index_in_S(Index_in_B);

% Reshaping the index by columns for idnetifying each reaction's transition
State_Transition_Index_Matrix = reshape(Index_in_S,nstates,number_reactions);
% fprintf('Transition matrix generated in %d seconds\n',toc)
end

% function State_Transition_Index_Matrix = State_Transition_Matrix_Finder(S,Stoichiometry)
% tic
% global nstates number_reactions
% State_Transition_Index_Matrix = zeros(nstates,number_reactions);
% for i = 1:number_reactions
%     % Here the Logic is - What are the new states from the current state
%     % due to different reactions. Another possible logic to code is -
%     % From which all state the system can reach current state.
%     
%     B = bsxfun(@plus,S',Stoichiometry(:,i))';
%     %% Removing unrealistic states
% %     B(any(B < 0,2),:) = nan;
%     %% ismember function is (very) slow which affects the speed of algorithm
%     [~,State_Transition_Index_Matrix(:,i)] = ismember(B,S,'rows');    
% end
% fprintf('Transition matrix generated in %d seconds\n',toc)
% end
