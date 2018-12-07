%% Model Name
CMEModel = sbiomodel('Huang1996 - Ultrasensitivity in MAPK cascade', 'Tag','Stochastic'); 

%% Species
addspecies(CMEModel,'MAPKKK activator (Ras)','InitialAmount',72265691','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPKKK inactivator','InitialAmount',722656903','InitialAmountUnits','molecule');
addspecies(CMEModel,'Mos','InitialAmount',7226569029','InitialAmountUnits','molecule');
addspecies(CMEModel,'Mos-P','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'Mek1','InitialAmount',2890627611360','InitialAmountUnits','molecule');
addspecies(CMEModel,'Mek1-P','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'Mek1-PP','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'Erk2','InitialAmount',2890627611360','InitialAmountUnits','molecule');
addspecies(CMEModel,'Erk2-P','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'Erk2-PP','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPK-Pase','InitialAmount',289062761136','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPKK-Pase','InitialAmount',722656903','InitialAmountUnits','molecule');
addspecies(CMEModel,'E1_Mos','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'E2_Mos-P','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'P-Mos_Mek1','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'P-Mos_P-Mek1','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'PP-Mek1_Erk2','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'PP-Mek1_P-Erk2','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPKK-Pase_PP-Mek1','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPKK-Pase_P-Mek1','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPK-Pase_PP-Erk2','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'MAPK-Pase_P-Erk2','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'K_PP_norm','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'KK_PP_norm','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'KKK_P_norm','InitialAmount',0','InitialAmountUnits','molecule');
addspecies(CMEModel,'relative maximal K_PP','InitialAmount',0','InitialAmountUnits','molecule');

%% Boundary Condition
BoundaryCondition = zeros(1,26);
% BoundaryCondition = [];

%% Reactions 
addreaction(CMEModel,'Mos + [MAPKKK activator (Ras)] -> E1_Mos','ReactionRate','4.1513e-10');
addreaction(CMEModel,'E1_Mos -> Mos + [MAPKKK activator (Ras)]','ReactionRate','150');
addreaction(CMEModel,'E1_Mos -> [MAPKKK activator (Ras)] + [Mos-P]','ReactionRate','150');
addreaction(CMEModel,'[Mos-P] + [MAPKKK inactivator] -> [E2_Mos-P]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[E2_Mos-P] -> [Mos-P] + [MAPKKK inactivator]','ReactionRate','150');
addreaction(CMEModel,'[E2_Mos-P] -> [MAPKKK inactivator] + Mos','ReactionRate','150');
addreaction(CMEModel,'Mek1 + [Mos-P] -> [P-Mos_Mek1]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[P-Mos_Mek1] -> Mek1 + [Mos-P]','ReactionRate','150');
addreaction(CMEModel,'[P-Mos_Mek1] -> [Mek1-P] + [Mos-P]','ReactionRate','150');
addreaction(CMEModel,'[Mek1-P] + [MAPKK-Pase] -> [MAPKK-Pase_P-Mek1]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[MAPKK-Pase_P-Mek1] -> [Mek1-P] + [MAPKK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[MAPKK-Pase_P-Mek1] -> Mek1 + [MAPKK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[Mek1-P] + [Mos-P] -> [P-Mos_P-Mek1]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[P-Mos_P-Mek1] -> [Mek1-P] + [Mos-P]','ReactionRate','150');
addreaction(CMEModel,'[P-Mos_P-Mek1] -> [Mek1-PP] + [Mos-P]','ReactionRate','150');
addreaction(CMEModel,'[Mek1-PP] + [MAPKK-Pase] -> [MAPKK-Pase_PP-Mek1]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[MAPKK-Pase_PP-Mek1] -> [Mek1-PP] + [MAPKK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[MAPKK-Pase_PP-Mek1] -> [Mek1-P] + [MAPKK-Pase]','ReactionRate','150');
addreaction(CMEModel,'Erk2 + [Mek1-PP] -> [PP-Mek1_Erk2]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[PP-Mek1_Erk2] -> Erk2 + [Mek1-PP]','ReactionRate','150');
addreaction(CMEModel,'[PP-Mek1_Erk2] -> [Erk2-P] + [Mek1-PP]','ReactionRate','150');
addreaction(CMEModel,'[Erk2-P] + [MAPK-Pase] -> [MAPK-Pase_P-Erk2]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[MAPK-Pase_P-Erk2] -> [Erk2-P] + [MAPK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[MAPK-Pase_P-Erk2] -> Erk2 + [MAPK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[Erk2-P] + [Mek1-PP] -> [PP-Mek1_P-Erk2]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[PP-Mek1_P-Erk2] -> [Erk2-P] + [Mek1-PP]','ReactionRate','150');
addreaction(CMEModel,'[PP-Mek1_P-Erk2] -> [Mek1-PP] + [Erk2-PP]','ReactionRate','150');
addreaction(CMEModel,'[Erk2-PP] + [MAPK-Pase] -> [MAPK-Pase_PP-Erk2]','ReactionRate','4.1513e-10');
addreaction(CMEModel,'[MAPK-Pase_PP-Erk2] -> [Erk2-PP] + [MAPK-Pase]','ReactionRate','150');
addreaction(CMEModel,'[MAPK-Pase_PP-Erk2] -> [Erk2-P] + [MAPK-Pase]','ReactionRate','150');
