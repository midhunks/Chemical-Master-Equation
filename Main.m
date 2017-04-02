%% Step 1: Add subdirectories to the search path
if ispc, b='\'; else, b='/'; end % defining forward/bckward slashes
addpath([pwd,b,'Generalized_Functions',b]);
addpath([pwd,b,'Models',b]);
addpath([pwd,b,'CME',b]);
addpath([pwd,b,'CME',b,'State_Space_Builder',b]);
addpath([pwd,b,'CME',b,'Dynamics_Builder',b]);

%% Step 2: Clearing variables and setting defaults for plots
Start_UP

%% Step 3: Generate the biochemical model
% File which describe the Biochemical model

MM

%% Step 4: Generate the CME

CME