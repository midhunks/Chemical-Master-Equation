%% Model Name
CMEModel = sbiomodel('Gene Toggle', 'Tag','Stochastic');

%% Species
addspecies(CMEModel, 'A', 'InitialAmount', 85e0,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'B', 'InitialAmount', 5e0,'InitialAmountUnits','molecule');

%% Reactions
addreaction(CMEModel, 'null -> A');
addreaction(CMEModel, 'A -> null');
addreaction(CMEModel, 'null -> B');
addreaction(CMEModel, 'B -> null');

%% Boundary condtion on molecular population of each species for open systems
% In an open system, there will be species with unbounded population.
% If there exists a conservation of some species, they will be aleady
% bounded. For unconserved species, the upper bound values should be
% defined so that a finite CME can be generated.

BoundaryCondition = [1e5,2e5,250];
