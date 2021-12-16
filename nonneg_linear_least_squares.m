%% Non Negative Least squares file 
%
% In this file we explore the use of non negative least squares to fit the
% one compartment model 


%% Non negative least squares method (using matlab)
% In this example we illustrate how we use the the inbuild matlab function
% for non negative least squares to estimate parameters to a single voxel 

[vY, vEchoTime] = GetVoxelMeasurements(Data, ET, 37, 35, 29); % Get the measurements and echo time values for a single value 
X = [ones(length(vY), 1) -vEchoTime]; % Set up the design matrix 
vParams = lsqnonneg(X, log(vY)); % Fit to the log of the signal 
fNNS0 = exp(vParams(1));         % Convert the log(S0) back to S0 
fNNT2 = 1/vParams(2);            % Convert the 1/T2 back to T2 
fResidual = OneCompSSD([fNNS0 fNNT2], vY, vEchoTime); % Compute the squared error value for our estimate
variance = fResidual/(length(X(:, 1)) - length(vParams));
fSTE = exp(X * vParams);
% [fS0s, fT2s, d] = BootStrap(vEchoTime, fSTE, variance);
% rangesfS0s = ConfidenceInterval(fS0s);
% rangesfT2s = ConfidenceInterval(fT2s);
%plot(vY, ' bs', 'MarkerSize', 16, 'LineWidth', 4);
%hold on;
%plot(fSTE, ' rx', 'MarkerSize', 16, 'LineWidth', 4);
% histogram(fT2s); 
% xlabel('T2');
% ylabel('Frequency');
%[fS0s, fT2s] = BootStrap(vEchoTime, fSTE, variance);
%% Linear Mapping of an entire brain 
% In this example we have encapsulated the method for non negative linear
% fitting into the NonnegLeastSquaresFit function. We then use this
% function to fit the entire brain. takes ~ 60secds

nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonNegS0Map = zeros(nDimX, nDimY, nDimZ); % Saves the S0 parameter mapping
mNonNegT2Map = zeros(nDimX, nDimY, nDimZ); % Saves the T2 parameter mapping
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 1:nDimZ
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice)) %Make sure we have a single mask that includes this voxel 
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = NonnegLeastSquaresFit(vY, vEchoTime);  
                mNonNegS0Map(nRow, nCol, nSlice) = vParams(1);
                mNonNegT2Map(nRow, nCol, nSlice) = vParams(2);
            end
        end
    end
end
toc;

%% Functions 
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
        vY = bootStrapData(:,i); 
        mX = [ones(nRows, 1) -vEchoTimes]; 
        vBeta = lsqnonneg(mX, log(vY));      
        fS0s(1, i) = exp(vBeta(1));
        fT2s(1, i) = (1/vBeta(2));
    end
end

