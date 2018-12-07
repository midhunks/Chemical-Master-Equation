%% Model Name
CMEModel = sbiomodel('LinearModel_3Species', 'Tag','Stochastic');

%% Species
addspecies(CMEModel, 'A', 'InitialAmount', 3e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'B', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'C', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');

%% Reactions
% addreaction(CMEModel, 'null -> A', 'ReactionRate', '1e0');
addreaction(CMEModel, 'A -> B', 'ReactionRate', '1e3');
addreaction(CMEModel, 'B -> A', 'ReactionRate', '5e2');
addreaction(CMEModel, 'B -> C', 'ReactionRate', '1e0');
addreaction(CMEModel, 'C -> B', 'ReactionRate', '0e0');

%% Boundary condtion on molecular population of each species for open systems
% In an open system, there will be species with unbounded population.
% If there exists a conservation of some species, they will be aleady
% bounded. For unconserved species, the upper bound values should be
% defined so that a finite CME can be generated.

BoundaryCondition = [];
