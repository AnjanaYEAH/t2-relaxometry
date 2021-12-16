%% Weighted Linear Least Squares 
%
% The aim of this script is to explore the results from the weighted linear
% least squares method. The key to this method is to compute an estimation
% for the noise of the measurements in the data. To do so, we get the error
% vector for an classic linear least squares fit, then proceed to compute
% the weight matrix. 
% Therefore this is a two part process where:
% 1. error vector using GLM error projection
% 2. using the error vector compute the weighted matrix 
% 3. Conduct the weighted least squares 


%% GLM to compute weighted least squares at a single voxel

% This section uses the standard GLM method with the project matrix to
% compute the error vector. Using this error vector we can compute the
% weight matrix. This enables us to compuate the weighted least squares
% parameter estimate.

% Design matrix 
[vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 37, 35, 29);
mX = [ones(length(vY), 1) vEchoTimes]; 
mPx = mX * pinv(mX' * mX) * mX';  
% Compute Rx = I - Px 
[nRows nCols] = size(mPx);
mRx = eye(nCols) - mPx;
% Compute the error estmation
vErrorEstimate = mRx * log(vY);
% Compute weight matrix 
mW = diag(1./vErrorEstimate);
% Compute parameter estimate using weighted least squares 
vParams = pinv(mX' * mW' * mW * mX) * (mX' * (mW'*mW)) * log(vY);  
fWS0 = exp(vParams(1));
fWT2 = -1/(vParams(2));
fWResidual = OneCompSSD([fWS0 fWT2], vY, vEchoTimes);

fSTE = exp(mX * vParams); 
variance = fWResidual/(length(mX(:, 1)) - length(vParams));
[fS0s, fT2s, data] = BootStrap(vEchoTimes, fSTE, variance);
rangesfS0s = ConfidenceInterval(fS0s);
rangesfT2s = ConfidenceInterval(fT2s);
histogram(fS0s); 
xlabel('S0');
ylabel('Frequency');

%% GLM to compute weighted least squares for the whole map
%
% In this section we encapsulate finding the parameters for a single voxel
% into a single function called LinearLeastSquares() takes ~ 90 seconds
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mS0Map = zeros(nDimX, nDimY, nDimZ);
mT2Map = zeros(nDimX, nDimY, nDimZ);
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 1:nDimZ
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = WeightedLeastSquaresFit(vY, vEchoTime); 
                mS0Map(nRow, nCol, nSlice) = vParams(1);
                mT2Map(nRow, nCol, nSlice) = vParams(2);
            end
        end
    end
end
toc;
disp(toc);






%% Functions 
function vParams = WeightedLeastSquaresFit(vY, vEchoTimes)
    mX = [ones(length(vY), 1) vEchoTimes]; 
    mPx = mX * pinv(mX' * mX) * mX';  
    % Compute Rx = I - Px 
    [nRows nCols] = size(mPx);
    mRx = eye(nCols) - mPx;
    % Compute the error estmation
    vErrorEstimate = mRx * log(vY);
    % Compute weight matrix 
    mW = diag(1./vErrorEstimate);
    % Compute parameter estimate using weighted least squares 
    vParams = pinv(mX' * mW' * mW * mX) * (mX' * (mW'*mW)) * log(vY);  
    fWS0 = exp(vParams(1));
    fWT2 = -1/(vParams(2));
    vParams = [fWS0 fWT2]; 
end 


function vParams = NonnegLeastSquaresFit(vY, vEchoTime)
    % Inputs:
    %      vY - n x 1 column vector of signal observations 
    %      vEchoTime - n x 1 column vector of corresponding echo time 
    % Output: 
    %      vParams - 2 x 1 column vector for the parameters S0 and T2
    % Uses the non negative Linear Least Squares method to estimate parameters 
    vNegativeIndices = find(vY <= 0); 
    while(length(vNegativeIndices) > 0)    %In this while loop we remove negative values from our data. 
        nIndex = vNegativeIndices(1); 
        vY(nIndex) = [];
        vEchoTime(nIndex) = []; 
        vNegativeIndices = find(vY <= 0);
    end 
    X = [ones(length(vY), 1) -vEchoTime];
    vParams = lsqnonneg(X, log(vY));
    fNNS0 = exp(vParams(1)); 
    fNNT2 = 1/vParams(2);
    vParams = [fNNS0 fNNT2]; 
end 

function vParams = LinearLeastSquaresFit(vY, vEchoTimes)
    % Inputs:
    %      vY - n x 1 column vector of signal observations 
    %      vEchoTime - n x 1 column vector of corresponding echo time 
    % Output: 
    %      vParams - 2 x 1 column vector for the parameters S0 and T2
    % Uses the basic Linear Least Squares method to estimate parameters 
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
    % Gets all the signal observations for a single voxel it their
    % corresponding echo time values. 
        % Find out how many observations there are for a single voxel for all
    % users 
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

function [fS0s, fT2s, bootStrapData] = BootStrap(vEchoTimes, S, variance)
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