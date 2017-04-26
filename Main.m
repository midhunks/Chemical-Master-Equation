%% Step 0: Add subdirectories to the search path
if ispc, b='\'; else, b='/'; end % defining forward/bckward slashes
addpath([pwd,b,'Generalized_Functions',b]);

%% Step 1: Clearing variables and setting defaults for plots
Start_UP

%% Step 2: Generate the biochemical model
% File which describe the Biochemical model

% MM
MAPK
% Two_Component_Pathway

%% Step 3: Generate the CME

CME