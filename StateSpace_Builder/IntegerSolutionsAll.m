function S = IntegerSolutionsAll(Moiety_Conservation, Conservation_Sum)
%% finding all combinations of linearly independent (n-k) species
    % The function identifies al possible states the system can reach
    % by means of conservation and upper bound.

%% Upper bound of molecular population of each species
% for a closed system, depending upon the conservation law and 
% initial molecular population, there is an upper bound on each species. 
% In open systems, we need to define the upper bound explicitely for
% generating a finite CME. For closed system algorithm  will automatically
% define the upper bopund of each species.

if ~isempty(Moiety_Conservation)    
    UB = floor(min(repmat(Conservation_Sum,1,size(Moiety_Conservation,2))./Moiety_Conservation,[],1));
else
    return
end

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
end