%% Model Name
CMEModel = sbiomodel('decaying-dimerizing reaction', 'Tag','Stochastic');

%% Species
addspecies(CMEModel, 'A', 'InitialAmount', 2e2,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'B', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'C', 'InitialAmount', 0e0,'InitialAmountUnits','molecule');

%% Reactions
% addreaction(CMEModel, 'null -> A', 'ReactionRate', '1e0');
addreaction(CMEModel, 'A -> null', 'ReactionRate', '1e0');
addreaction(CMEModel, '2 A -> B', 'ReactionRate', '2e-3');
addreaction(CMEModel, 'B -> 2 A', 'ReactionRate', '5e-1');
addreaction(CMEModel, 'B -> C', 'ReactionRate', '4e-2');

%% Boundary condtion on molecular population of each species for open systems
% In an open system, there will be species with unbounded population.
% If there exists a conservation of some species, they will be aleady
% bounded. For unconserved species, the upper bound values should be
% defined so that a finite CME can be generated.

BoundaryCondition = [1e5,2e5,250];
