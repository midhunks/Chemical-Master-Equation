%% Home folder
if ~exist('CME_Folder','var')
    % Use one of the following two lines
    CME_Folder = uigetdir('','Choose CME Folder');
    % CME_Folder = '/Chemical-Master-Equation';  needs to give accurate path
end

% Update CME
cd(CME_Folder);
system('git pull https://github.com/midhunks/Chemical-Master-Equation.git');

%% Step 1: Model Setup
% Run the CMEModel from Matlab (.m) model file. Template files are
% available in 'Models' folder.
cd(CME_Folder);  cd('Models');
[FileName,path] = uigetfile({'*.m','Model Files (*.m)';
                    '*.*',  'All Files (*.*)'}, 'Choose the model file');
run(FileName)

% Or simply use the following line using the name of the model file
% run('MM.m')

%% %% Step 2: Generate State space of the CME
cd(CME_Folder);     cd('StateSpace_Builder');
State_Space = StateSpace_Builder(CMEModel, BoundaryCondition, 0);

%% %% Step 3: Generatethe transition matrix of the CME
cd(CME_Folder);     cd('TransitionMatrix_Builder');
[Transition_Matrix, State_Transition_Index_Matrix, SSA_propensity_matrix]...
                            = TransitionMatrix_Builder(CMEModel,State_Space);
