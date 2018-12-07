% Model Name
CMEModel = sbiomodel('Schlogl', 'Tag','Stochastic');

% Species
addspecies(CMEModel, 'A', 'InitialAmount', 1e5,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'B', 'InitialAmount', 2e5,'InitialAmountUnits','molecule');
addspecies(CMEModel, 'X', 'InitialAmount', 25e1,'InitialAmountUnits','molecule');

% Reactions
addreaction(CMEModel, 'null -> A', 'ReactionRate', '1e0');
addreaction(CMEModel, 'A + 2 X  -> 3 X', 'ReactionRate', '3e-7');
addreaction(CMEModel, '3 X -> A + 2 X', 'ReactionRate', '1e-4');
addreaction(CMEModel, 'B -> X', 'ReactionRate', '1e-3');
addreaction(CMEModel, 'X -> B', 'ReactionRate', '3.5e0');

% Boundary condtion on molecular population of each species for open systems
In an open system, there will be species with unbounded population.
If there exists a conservation of some species, they will be aleady
bounded. For unconserved species, the upper bound values should be
defined so that a finite CME can be generated.

BoundaryCondition = [1e5,2e5,250];
