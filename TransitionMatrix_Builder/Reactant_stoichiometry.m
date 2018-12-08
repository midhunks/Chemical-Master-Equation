function S = Reactant_stoichiometry(Model)
%% Reactant Stoichimetry: By setting stoichiometry in terms of reactants
% In some reactions some species present as reactant and product.
% For e.g., A + B -> A + C. where A is a reactant and a product which act as a catalyst.
% Net gain of A is zero molecules. This will remove A from original
% Stoichiometry matrix and cannot idnetify the reactants properly which is
% necassary for computing reaction propensities in CME
ModelTemp = copyobj(Model);
number_reactions = length(ModelTemp.Reactions);
for i = 1:number_reactions
    Products = ModelTemp.Reactions(i).Products;
    if ~isempty(Products)
        rmproduct(ModelTemp.Reactions(i),Products);
    end
end
S = full(getstoichmatrix(ModelTemp));
end
