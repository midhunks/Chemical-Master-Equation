%% Setting global variables. (Please do not touch this secton !!!)
global minimal_CME Reaction_rates Conservation_temp ...
       Initial_Molecular_Population Boundary_condition
   
%% Minimal CME generation (Choose 0 or 1)
minimal_CME = 0;

%% Stoichiometry (Input as per the system) (Species order - R,L, RL, P, P*)
%   Model from  Brian Ingalls' notes Section 6.1.1

Stoichiometry=[-1  1  0  0;
               -1  1  0  0;
                1 -1  0  0;
                0  0 -1  1;
                0  0  1 -1];

%% Identify reactants: By setting stoichiometry in terms of reactants 
%	In some reactions some species present as reactant and product. 
%   For e.g., A + B -> A + C. where A is a reactant and a product. 
%   Net gain of A is zero molecules. This will remove A from original 
%   Stoichiometry matrix and cannot idnetify the reactants properly 
%   which is necassary for computing reaction propensities in CME.

Reactants_stoichiometry=[-1  0  0  0;
                         -1  0  0  0;
                          0 -1 -1  0;
                          0  0 -1  0;
                          0  0  0 -1];     

%% Reaction rates (equations or constants)
Reaction_rates =[5e0 1e0 1e3 3e3]; %Parameters of reaction_rate

%   The following function will generate an array of 'Reaction_rates'.
%   However, when the 'Reaction_rates' are changing with respect to the
%   current population of species, edit the following function accordingly.
%   Leave it commented if the rates are constants as before.

% global Rate
% Rate = Reaction_rate_Function(Reaction_rates,S);

%% Initial molecular population
% volume=1e-19;
% Avogadro=6.022140857e23;
% C_1=2e0;
% N_1=round(C_1*Avogadro*volume);
% C_4=8e0;
% N_4=round(C_4*Avogadro*volume);
% Initial_Molecular_Population=[N_1 1 0 N_4 0];
Initial_Molecular_Population=10*[1 1 1 1 1];
%% Boundary condtion on molecular population of each species for open systems
%   In an open system, there will be species with unbounded population.
%   If there exists a conservation of some species, they will be aleady
%   bounded. For unconserved species, the upper bound values should be
%   defined so that a finite CME can be generated. (Default) boundary
%   condition is set to 10 molecules each for (unconserved) species.
%   Conserved species will accept values irrespective of the input

Boundary_condition=[inf inf inf inf inf]; 

%% Conservation Laws
%   Keep this untouched as long as the conservation law identified from rref has
%   integer values. Otherwise, modify manually with integer conservation
%   values. First identify conservation law by running (uncomment the following line)

% Conservation = Conservation_Law (Stoichiometry)

%   Scale the result so that all columns have integer values.

Conservation_temp=[];
