%% Linear Least Squares fitting 
% 
% In this file we explore using the basic Linear least Squares fitting
% method to estimate our parameters 

%% Parameter modelling Linear Fit 
% In this section we get the measurements and echo times for voxel (30, 30,
% 30) and use it to estimate the parameters S0 and T2
% 
% [vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 30, 30, 30); %Get the observations and measurements 
% [nRows, nCols] = size(vY);
% vLogY = log(vY);  % Get the log of the observations 
% mX = [ones(nRows, 1) vEchoTimes]; %design matrix 
% vBeta = pinv(mX) * vLogY;         %parameter estimation using pinv 
% fS0 = exp(vBeta(1));              %converting log(S0) back to S0
% fT2 = -(1/vBeta(2));              %converting -1/T2 back to T2
% fResidual = OneCompSSD([fS0 fT2], vY, vEchoTimes); % Compute the squared error value for our estimated parameters
% 
% 
% %CSF region
% [vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 30, 55, 49);
% [nRows, nCols] = size(vY);
% vLogY = log(vY); 
% mX = [ones(nRows, 1) vEchoTimes];  
% vBeta = pinv(mX) * vLogY;         
% fS0 = exp(vBeta(1));            
% fT2 = -(1/vBeta(2));              
% fResidual = OneCompSSD([fS0 fT2], vY, vEchoTimes);
% fSTE = exp(mX * vBeta); 
% variance = fResidual/(length(mX(:, 1)) - length(vBeta));
% [fS0s, fT2s] = BootStrap(vEchoTimes, fSTE, variance);
% rangesfS0s = ConfidenceInterval(fS0s);
% rangesfT2s = ConfidenceInterval(fT2s);

%WM region 
[vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 35, 37, 29);
[nRows, nCols] = size(vY);
vLogY = log(vY); 
mX = [ones(nRows, 1) vEchoTimes];  
vBeta = pinv(mX) * vLogY;         
fS0 = exp(vBeta(1));            
fT2 = -(1/vBeta(2));              
fResidual = OneCompSSD([fS0 fT2], vY, vEchoTimes);

% fSTE = exp(mX * vBeta); 
% variance = fResidual/(length(mX(:, 1)) - length(vBeta));
% [fS0s, fT2s] = BootStrap(vEchoTimes, fSTE, variance);
% rangesfS0s = ConfidenceInterval(fS0s);
% rangesfT2s = ConfidenceInterval(fT2s);
% histogram(fT2s); 
% xlabel('T2');
% ylabel('Frequency');
%% Linear Mapping of an entire brain 
%
% In this section we encapsulate finding the parameters for a single voxel
% into a single function called LinearLeastSquares() takes ~ 50 seconds
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mS0Map = zeros(nDimX, nDimY, nDimZ);
mT2Map = zeros(nDimX, nDimY, nDimZ);
tic;
for nRow = 1:nDimX
    for nCol = 1:nDimX
        for nSlice = 1:nDimZ
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = LinearLeastSquaresFit(vY, vEchoTime); 
                mS0Map(nRow, nCol, nSlice) = vParams(1);
                mT2Map(nRow, nCol, nSlice) = vParams(2);
            end
        end
    end
end
toc;


%% Functions 
function vParams = LinearLeastSquaresFit(vY, vEchoTimes)
    % Inputs:
    %      vY - n x 1 column vector of signal observations 
    %      vEchoTime - n x 1 column vector of corresponding echo time 
    % Output: 
    %      vParams - 2 x 1 column vector for the parameters S0 and T2
    
    vNegativeIndices = find(vY <= 0); 
    while(length(vNegativeIndices) > 0)
        nIndex = vNegativeIndices(1); 
        vY(nIndex) = [];
        vEchoTimes(nIndex) = []; 
        vNegativeIndices = find(vY <= 0);
    end 
    [nRows, nCols] = size(vY); 
    vLogY = log(vY); 
    mX = [ones(nRows, 1) vEchoTimes]; %design matrix 
    vBeta = pinv(mX) * vLogY; 
    fS0 = exp(vBeta(1));
    fT2 = -(1/vBeta(2)); 
    vParams = [fS0 fT2]; 
end

function [vY, vEchoTime]  = GetVoxelMeasurements(Data, EchoTimes, nX, nY, nZ)
    % Input:
    %     "Data" cell array - The collection of all the 4D qt2 matrices
    %     "ET" cell array - The collection of all the Echo Time data 
    %      nX, nY, nZ - x, y and z coordinate of the voxel of interest
    % Output:
    %      vY - n x 1 column vector of all the signal observations from all
    %           6 cases at a specific voxel 
    % Gets all the signal observations for a single voxel and their
    % corresponding echo time values. 
    
    [dim1, dim2] = size(Data); 
    if(dim1 > dim2)
        nSubjects = dim1;
    else 
        nSubjects = dim2; 
    end 
    
    vEchoTime = cell2mat(EchoTimes(1)); 
    mSubjectData01 = cell2mat(Data(1)); 
    vY = squeeze(mSubjectData01(nX, nY, nZ, :)); 
    for s = 2:nSubjects
        vSubjectEchoTimes = cell2mat(EchoTimes(s)); 
        vEchoTime = vertcat(vEchoTime, vSubjectEchoTimes);
        %Now we collect all the obersvations 
        mSubjectData = cell2mat(Data(s)); 
        vY = vertcat(vY, squeeze(mSubjectData(nX, nY, nZ, :)));
        vY = double(vY);
    end 
end

function bMasked = GetGlobalMaskVal(Mask, nX, nY, nZ) 
    % Input 
    %      "Mask" - a 4D matrix that holds every mask 
    %      nX, nY, nZ - The x, y and z coordinates of the 
    % Output
    %      boolean bMasked - represents whether one of the masks includes
    %                        the voxel of interest
    
    nGlobalMaskVal = sum(Mask(nX, nY, nZ, :), 'all'); 
    if(nGlobalMaskVal > 0) 
        bMasked = 1; 
    else 
        bMasked = 0; 
    end 
end 

function [vY, vEchoTime] = GetSingleSubjectVoxelMeasurements(Data, EchoTimes, nSubject, nX, nY, nZ)
    % Works the same as GetVoxelMeasurements(Data, EchoTimes, nX, nY, nZ)
    % but gets the voxel signal observations for a single subject instead
    % of all the subjects. Hence in the 3 parameter we specify which
    % subject we want. 
    % 
    % Input:
    %     "Data" cell array - The collection of all the 4D qt2 matrices
    %     "ET" cell array - The collection of all the Echo Time data 
    %      nSubject (int) - The subject of interest
    %      nX, nY, nZ - x, y and z coordinate of the voxel of interest
    % Output:
    %      vY - n x 1 column vector of all the signal observations from all
    %           6 cases at a specific voxel 
    % Gets all the signal observations for a single voxel and their
    % corresponding echo time values. 
    vEchoTime = cell2mat(EchoTimes(nSubject)); 
    mSubjectData = cell2mat(Data(nSubject)); 
    vY = squeeze(mSubjectData(nX, nY, nZ, :)); 
    vY = double(vY);
end 

function vEstY = GetSyntheticSignal(vParams, vEchoTime)
    % Computes the Synthetic signal using our estimated parameters S0 and
    % T2 given the echo time 
    %
    % Input 
    %      vParams - 2 x 1 column vector that contains an estimated S0 and
    %                T2.
    %      vEchoTime - n x 1 column vector of echo times 
    % Output 
    %      vEstY - n x 1 column vector of synthetic signal 
    fS0 = vParams(1); 
    fT2 = vParams(2); 
    [nRows, nCols] = size(vEchoTime); 
    vEstY = zeros(nRows, 1); 
    for nRow = 1:nRows
        vEstY(nRow) = fS0 * exp(vEchoTime(nRow) * (-1.0 / fT2));
    end
    vEstY = double(vEstY);
end 

function [fS0s, fT2s, data] = BootStrap(vEchoTimes, S, variance)
    [nRows, nCols] = size(S);
    bootStrapSize = 1000;
    bootStrapData = zeros(length(S), bootStrapSize);
    for i = 1: length(S)
        for j = 1: bootStrapSize
            bootStrapData(i, j) = S(i) + normrnd(0, sqrt(200));
        end
    end
    fS0s = zeros(1, bootStrapSize);
    fT2s = zeros(1, bootStrapSize);
    for i = 1: bootStrapSize
        vLogY = log(bootStrapData(:,i)); 
        mX = [ones(nRows, 1) vEchoTimes]; 
        vBeta = pinv(mX) * vLogY;          
        fS0s(1, i) = exp(vBeta(1));
        fT2s(1, i) = -(1/vBeta(2));
    end
end