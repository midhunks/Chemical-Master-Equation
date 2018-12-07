%% Model Name
CMEModel = sbiomodel('Michaelis-Menten', 'Tag','Stochastic');

%% Species
addspecies(CMEModel, 'S', 'InitialAmount', 3e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'E', 'InitialAmount', 3e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'C', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'P', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');

%% Reactions
% addreaction(CMEModel, 'null -> S', 'ReactionRate', '1e0');
addreaction(CMEModel, 'S + E -> C', 'ReactionRate', '1e3');
addreaction(CMEModel, 'C -> E + S', 'ReactionRate', '5e2');
addreaction(CMEModel, 'C -> E + P', 'ReactionRate', '5e0');
addreaction(CMEModel, 'E + P -> C', 'ReactionRate', '1e0');

%% Boundary condtion on molecular population of each species for open systems
% In an open system, there will be species with unbounded population.
% If there exists a conservation of some species, they will be aleady
% bounded. For unconserved species, the upper bound values should be
% defined so that a finite CME can be generated.

x = 4;
BoundaryCondition = [x, 3, 3, x];
