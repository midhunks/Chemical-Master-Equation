function [Index] = Combvec_Optimization_index(Moeity_Conservation,UB)
%% While using combvec, we have to keep the vectors generated using combvec
    % less exploding. In every iteration, we will generate a set 
    % of column vectors using combvec and then remove some rows which are 
    % not satisfying the Moeity_Conservation laws. If the first vectors are large 
    % enough, then by iteratively removing the rows using Moeity_Conservation, 
    % the newer vectors formed with olde combvec output will be less expoloding. 
    % The following re-indexing will help it.

% global number_species 
% dim_cons = size(Moeity_Conservation,2);
% Index = 1:number_species;
% Logical_Moeity_Conservation = (Moeity_Conservation ~= 0); %Generate matrix with only ones and zeros
% 
% % Re-ndexing as per the descending order of the upper bound
% [~,m] = sort(UB,'descend');
% Index = Index(m);
% Logical_Moeity_Conservation = Logical_Moeity_Conservation(m,:);
% 
% % Identify and remove the dependent sepcies
% Dependent_Species_Index = find(sum(Logical_Moeity_Conservation,2)~=1);
% Logical_Moeity_Conservation(Dependent_Species_Index,:)= [];
% Index(Dependent_Species_Index)=[];
% 
% % Identifying the independent rows with largest upperbound
% for i = 1:dim_cons
%     Temp = Logical_Moeity_Conservation(i,:);    
%     Temp_Repetition_index = ismember(Logical_Moeity_Conservation, Temp,'rows');
%     Temp_Repetition_index(i) = 0; % To avoid the removal of parent Temp
% 
%     Logical_Moeity_Conservation(Temp_Repetition_index,:)= [];
%     Index(Temp_Repetition_index) = [];
% end
% 
% Index = [Index setdiff(1:number_species,Index)];


%%
Species_Index = 1:size(Moeity_Conservation,1)

% [~,m] = sort(UB,'descend')
% Moeity_Conservation = Moeity_Conservation(m,:)

[Moeity_Conservation_Temp, Independent_Species_Index] = rref(Moeity_Conservation')

[~,Index] = sort(UB(Independent_Species_Index),'descend')
Independent_Species_Index = Independent_Species_Index(Index)
Index = [Independent_Species_Index setdiff(Species_Index,Independent_Species_Index)]
end




