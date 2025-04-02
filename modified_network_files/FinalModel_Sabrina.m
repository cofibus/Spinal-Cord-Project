%% FULL NETWORK WITHOUT GENETIC BIAS (SPATIAL BIASES)
clear
MyNetwork = NetworkParameters();
g = 2;
gc = 0.5;

% Ventral Interneurons
MyNetwork.AddCelltype(Celltype.V2a_2,'bi','ven','lat','','ipsi','','',g*3500,'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V2a_1,'cau','ven','','','ipsi','','',500,'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V1_Foxp2,'ro','ven','lat','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Ren,'loc','ven','loc','','ipsi','','',g*500,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Pou6f2,'ro','loc','lat','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Sp8,'ro','loc','lat','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V2b,'cau','ven','lat','','ipsi','','',g*2000,'SynStrength',gc*1,'Gain',1)
MyNetwork.AddCelltype(Celltype.V0d,'loc','ven','med','','contra','','',g*250,'SynStrength',gc*1,'Gain',1) 
MyNetwork.AddCelltype(Celltype.V0v,'ro','loc','med','','contra','','',g*1700,'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.DI6,'cau','dor','med','','contra','','',g*2500,'SynStrength',gc*1,'Gain',0.1) 
MyNetwork.AddCelltype(Celltype.V3,'cau','ven','lat','','contra','','',500,'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.MN,'loc','ven','lat','','ipsi','','',300,'SynStrength',gc*1,'Gain',1);

MyNetwork.SetCellBiasPop('MN','All',0.1);
N = Network([append('T',string([13])),append('L',string([1:6])),append('S',string([1:4]))],1,1,MyNetwork,'UseGeneticBias',0);
N_wo = N;
save('N_wo.mat', 'N_wo', '-v7.3');
%% FULL NETWORK WITH GENETIC BIAS
clear
MyNetwork = NetworkParameters();
g = 2;
gc = 0.5;
% Ventral Interneurons
MyNetwork.AddCelltype(Celltype.V2a_2,'bi','','','','ipsi','','',g*3500,'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V2a_1,'cau','','','','ipsi','','',500,'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V1_Foxp2,'ro','','','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Ren,'loc','','','','ipsi','','',g*500,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Pou6f2,'ro','','','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Sp8,'ro','','','','ipsi','','',g*1000,'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V2b,'cau','','','','ipsi','','',g*2000,'SynStrength',gc*1,'Gain',1)
MyNetwork.AddCelltype(Celltype.V0d,'loc','','','','contra','','',g*250,'SynStrength',gc*1,'Gain',1) 
MyNetwork.AddCelltype(Celltype.V0v,'ro','','','','contra','','',g*1700,'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.DI6,'cau','','','','contra','','',g*2500,'SynStrength',gc*1,'Gain',0.1) 
MyNetwork.AddCelltype(Celltype.V3,'cau','','','','contra','','',500,'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.MN,'loc','','','','ipsi','','',300,'SynStrength',gc*1,'Gain',1);

MyNetwork.SetCellBiasPop('MN','All',0.1);
N = Network([append('T',string([13])),append('L',string([1:6])),append('S',string([1:4]))],1,1,MyNetwork,'UseGeneticBias',1);
N_w = N;
save('N_w.mat', 'N_w', '-v7.3');

%% FULL NETWORK RANDOM PROPERLY RANDOMIZE!!!!
clear
MyNetwork = NetworkParameters();
g = 2;
gc = 0.5;

% Ventral Interneurons
MyNetwork.AddCelltype(Celltype.V2a_2,'bi','ven','lat','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V2a_1,'cau','ven','','','ipsi','','',randi(3000),'SynStrength',gc*1,'Gain',1);
MyNetwork.AddCelltype(Celltype.V1_Foxp2,'ro','ven','lat','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Ren,'loc','ven','loc','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Pou6f2,'ro','loc','lat','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V1_Sp8,'ro','loc','lat','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1)
MyNetwork.AddCelltype(Celltype.V2b,'cau','ven','lat','','ipsi','','',g*randi(3000),'SynStrength',gc*1,'Gain',1)
MyNetwork.AddCelltype(Celltype.V0d,'loc','ven','med','','contra','','',g*randi(3000),'SynStrength',gc*1,'Gain',1) 
MyNetwork.AddCelltype(Celltype.V0v,'ro','loc','med','','contra','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.DI6,'cau','dor','med','','contra','','',g*randi(3000),'SynStrength',gc*1,'Gain',0.1) 
MyNetwork.AddCelltype(Celltype.V3,'cau','ven','lat','','contra','','',randi(3000),'SynStrength',gc*1,'Gain',0.1); 
MyNetwork.AddCelltype(Celltype.MN,'loc','ven','lat','','ipsi','','',randi(3000),'SynStrength',gc*1,'Gain',1);

MyNetwork.SetCellBiasPop('MN','All',0.1);
N = Network([append('T',string([13])),append('L',string([1:6])),append('S',string([1:4]))],1,1,MyNetwork,'UseGeneticBias',0);
N_r = N;
save('N_r.mat', 'N_r', '-v7.3');
%% RUN SIMULATION
t_steps = 10000;
[~,IxT] = ismember(N.Types,N.Parameters.Types);
I_e = zeros(size(N.ConnMat,1),t_steps); 
I_e(:) = 20;
N.Simulate(t_steps,'I_e',I_e); 