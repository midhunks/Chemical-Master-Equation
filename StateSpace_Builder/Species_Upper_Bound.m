function UB = Species_Upper_Bound(Moeity_Conservation, Initial_Molecular_Population)
%% Upper bound of molecular population of each species
% for a closed system, depending upon the conservation law and
% initial molecular population, there is an upper bound on each species.
% In open systems, we need to define the upper bound explicitely for
% generating a finite CME. For closed system algorithm  will automatically
% define the upper bopund of each species.

% Upper bound of ith species from the minimum of conservation sum
% where ith species is present.
if ~isempty(Moeity_Conservation)
    Conservation_Sum = Initial_Molecular_Population*Moeity_Conservation;
    UB = min(floor(repmat(Conservation_Sum,size(Moeity_Conservation,1),1)./Moeity_Conservation)');
end

% global number_species
% UB = zeros(1,number_species);
% % checks whether there is any conservation law exists
% if isempty(Moeity_Conservation)
%     if isempty(Boundary_condition)
%         warning on
%         warning('Open species exists. Define upper bound for open-species.')
%         warning('You can define Boundary condition while the model creation time.')
%         warning('Usage: Add the following lines in the model input file')
%         fprintf('global Boundary_condition;\n')
%         fprintf('Boundary_condition = zeros(number_species,1);\n')
%         fprintf(strcat('Boundary_condition([',num2str(1:number_species),']) = ["Upper bound of open species here"];\n'))
%     end
%     %upper bound of all species is set from the boundary condition
%     UB = Boundary_condition;
% else
%     % Closed (Conserved) Species' upper bound
%     Closed_Species_index = find(sum(Moeity_Conservation,2) ~= 0);
%     % Upper bound of ith species from the minimum of conservation sum
%     % where ith species is present.    
%     if ~isempty(Closed_Species_index)
%         Conservation_temp = Moeity_Conservation(Closed_Species_index,:);
%         Conservation_Sum = Initial_Molecular_Population(Closed_Species_index)*Conservation_temp;
%         UB(Closed_Species_index) = min(floor(repmat(Conservation_Sum,length(Closed_Species_index),1)./Conservation_temp)');
%     end
%         
%     % Open (Unconserved) Species' upper bound
%     Open_Species_index = find(sum(Moeity_Conservation,2) == 0);
%     if ~isempty(Open_Species_index)
%         %upper bound of ith species is set from the boundary condition
%         if isempty(Boundary_condition)
%             warning on
%             warning('Open species exists. Define upper bound for open-species.')
%             warning('You can define Boundary condition while the model creation time.')
%             warning('Usage:Add the following lines in the model input file')
%             fprintf('global Boundary_condition;\n')
%             fprintf('Boundary_condition = zeros(number_species,1);\n')
%             fprintf(strcat('Boundary_condition([',num2str(Open_Species_index),']) = ["Upper bound of open species here"];\n'))
%             
%             Boundary_condition = zeros(number_species,1);
%             for i = 1:length(Open_Species_index)
%                 prompt = strcat('Upper bound of species- ', num2str(Open_Species_index(i)),' :  ');
%                 Boundary_condition(Open_Species_index(i)) = input(prompt);
%             end
%         end
%         UB(Open_Species_index) = Boundary_condition(Open_Species_index);
%     end
% end
%%
end