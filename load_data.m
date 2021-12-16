%% Main data
% Load in all the Nifti datasets 

mCase01Mask = niftiread("./project3_data0/case01-mask.nii.gz");
mCase01qt2 = niftiread("./project3_data0/case01-qt2_reg.nii.gz");
mCase01seg = niftiread("./project3_data0/case01-seg.nii.gz");

mCase02Mask = niftiread("./project3_data0/case02-mask.nii.gz");
mCase02qt2 = niftiread("./project3_data0/case02-qt2_reg.nii.gz");
mCase02seg = niftiread("./project3_data0/case01-seg.nii.gz");

mCase03Mask = niftiread("./project3_data0/case03-mask.nii.gz");
mCase03qt2 = niftiread("./project3_data0/case03-qt2_reg.nii.gz"); 
mCase03seg = niftiread("./project3_data0/case01-seg.nii.gz");

mCase04Mask = niftiread("./project3_data1/case04-mask.nii.gz");
mCase04qt2 = niftiread("./project3_data1/case04-qt2_reg.nii.gz");
mCase04seg = niftiread("./project3_data0/case01-seg.nii.gz");

mCase05Mask = niftiread("./project3_data1/case05-mask.nii.gz");
mCase05qt2 = niftiread("./project3_data1/case05-qt2_reg.nii.gz"); 
mCase05seg = niftiread("./project3_data0/case01-seg.nii.gz");

mCase06Mask = niftiread("./project3_data1/case06-mask.nii.gz");
mCase06qt2 = niftiread("./project3_data1/case06-qt2_reg.nii.gz"); 
mCase06seg = niftiread("./project3_data0/case01-seg.nii.gz");

%% Echo time vectors
%Extract into echo times 
filePtr = fopen("./project3_data0/case01-TEs.txt"); 
vCase01ET = textscan(filePtr, "%f");
vCase01ET = cell2mat(vCase01ET);

filePtr = fopen("./project3_data0/case02-TEs.txt"); 
vCase02ET = textscan(filePtr, "%f"); 
vCase02ET = cell2mat(vCase02ET); 

filePtr = fopen("./project3_data0/case03-TEs.txt"); 
vCase03ET = textscan(filePtr, "%f"); 
vCase03ET = cell2mat(vCase03ET); 

filePtr = fopen("./project3_data1/case04-TEs.txt"); 
vCase04ET = textscan(filePtr, "%f"); 
vCase04ET = cell2mat(vCase04ET); 

filePtr = fopen("./project3_data1/case05-TEs.txt"); 
vCase05ET = textscan(filePtr, "%f"); 
vCase05ET = cell2mat(vCase05ET); 

filePtr = fopen("./project3_data1/case06-TEs.txt"); 
vCase06ET = textscan(filePtr, "%f"); 
vCase06ET = cell2mat(vCase06ET); 


%% Data cell structure that group up all 6 users 
Data{1} = mCase01qt2; 
Data{2} = mCase02qt2;
Data{3} = mCase03qt2; 
Data{4} = mCase04qt2; 
Data{5} = mCase05qt2; 
Data{6} = mCase06qt2; 

% Segmentation data structure that group up all 6 users
Dataseg{1} = mCase01seg;
Dataseg{2} = mCase02seg;
Dataseg{3} = mCase03seg;
Dataseg{4} = mCase04seg;
Dataseg{5} = mCase05seg;
Dataseg{6} = mCase06seg;

% Stores all the masks into a 4D matrix 
Masks = zeros(96, 96, 55, 6); 
Masks(:, :, :, 1) = mCase01Mask; 
Masks(:, :, :, 2) = mCase02Mask; 
Masks(:, :, :, 3) = mCase03Mask; 
Masks(:, :, :, 4) = mCase04Mask; 
Masks(:, :, :, 5) = mCase05Mask; 
Masks(:, :, :, 6) = mCase06Mask; 

% Store all the Echo times into cell arrays 
ET{1} = vCase01ET; 
ET{2} = vCase02ET; 
ET{3} = vCase03ET; 
ET{4} = vCase04ET; 
ET{5} = vCase05ET; 
ET{6} = vCase06ET; 





