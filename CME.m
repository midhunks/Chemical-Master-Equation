%% Home folder
CME_Folder = uigetdir('','Choose CME Folder');
% OR you can complete path if you are using the following line
% CME_Folder = '/Chemical-Master-Equation'; 
cd(CME_Folder);

%% Step 1: Model Setup
% Gnerate the Matlab model file from SBML (.xml) model file
cd(CME_Folder);  cd('Models');

% Run the CMEModel from Matlab (.m) model file
[FileName,path] = uigetfile({'*.m','Model Files (*.m)';
    '*.*',  'All Files (*.*)'}, ...
    'Choose the model file');
run(FileName)

% Or simply use the following line
% run('MM.m')

%% %% Step 2: Generate State space of the CME
cd(CME_Folder);     cd('State_Space_Builder');
State_Space = State_Builder(Stoichiometry);

%% %% Step 3: Generatethe transition matrix of the CME
cd(CME_Folder);     cd('TransitionMatrix_Builder');
[Transition_Matrix, State_Transition_Index_Matrix,SSA_propensity_matrix]...
        = Dynamics_Builder(State_Space, Stoichiometry, Reactants_stoichiometry);
