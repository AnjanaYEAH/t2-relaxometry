%% Non linear fitting of one compartment 
% We will use the non negative least squares as initial staring point for
% our iterative approach of finding the minimum 


%Define various options for the non-linear fitting algorithm 
h = optimset('MaxFunEvals', 2000, 'Algorithm', 'quasi-newton','Display', 'off', 'TolX', 1e-6,'TolFun', 1e-6);


%% First I need to find out how many times I need to iterate before I get the global minimum 
[vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 57, 12, 29);
startx = [mNonNegS0Map(57, 12, 29) mNonNegT2Map(57, 12, 29)]; 
[vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fminunc('OneCompSSD', startx, h, vY, vEchoTimes);
fS0 = vNonLinearParams(1);
fT2 = vNonLinearParams(2);
fSTE = fS0 * exp(-vEchoTimes / fT2);
% variance = RESNORM/(length(vEchoTimes) - length(startx));
% [fS0s, fT2s, d] = BootStrap(vEchoTimes, fSTE, startx, variance, h);
% rangesfS0s = ConfidenceInterval(fS0s);
% rangesfT2s = ConfidenceInterval(fT2s);


nIter = 100; 
vResnorms = zeros(nIter, 1); 
vParams = zeros(nIter, 2); 
vResnorms(100) = RESNORM;
vParams(100, :) = vNonLinearParams; 
for i = 1:nIter-1 
    fNewS0Start = mNonNegS0Map(57, 12, 29) + randn()*500; % linear square fitting S0 = ~ 1000s
    fNewT2Start = mNonNegT2Map(57, 12, 29) + randn()*50;  % linear square fitting T2 = ~ 70s 
    startx = [fNewS0Start fNewT2Start]; 
    [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fminunc('OneCompSSD', startx, h, vY, vEchoTimes);
    disp(RESNORM);
    vResnorms(i) = RESNORM; 
    vParams(i, :) = vNonLinearParams; 
end

%% Find the number of global minimums found
vRoundedResnorm = round(vResnorms, 4); 
nNumGlobalMins = length(find(vRoundedResnorm == min(vRoundedResnorm))); 
disp(nNumGlobalMins);

%% Fmincon 
% [vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 37, 35, 29);
% startx = [mNonNegS0Map(37, 35, 29); mNonNegT2Map(37, 35, 29)]; 
% % Ax <= b 
% % x = [S0 T2] 
% % A*x = [-S0; -T2]  <= b = [0; 0] 
% A = [-1 0; 0 -1];
% b = [0; 0];
% f = @(startx)OneCompSSD(startx, vY, vEchoTimes);
% [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
% variance = RESNORM/(length(vEchoTimes) - length(startx));
% fS0 = vNonLinearParams(1);
% fT2 = vNonLinearParams(2);
% fSTE = fS0 * exp(-vEchoTimes / fT2);
% [fS0s, fT2s, d] = BootStrap2(vEchoTimes, fSTE, startx, variance);
% rangesfS0s = ConfidenceInterval(fS0s);
% rangesfT2s = ConfidenceInterval(fT2s);
% histogram(fT2s); 
% xlabel('T2');
% ylabel('Frequency');
% 
% nIter = 100; 
% vResnorms = zeros(nIter, 1); 
% vParams = zeros(nIter, 2); 
% vResnorms(100) = RESNORM;
% vParams(100, :) = vNonLinearParams;
% for i = 1:nIter-1 
%     fNewS0Start = mNonNegS0Map(30, 30, 30) + randn()*500; 
%     fNewT2Start = mNonNegT2Map(30, 30, 30) + randn()*500;
%     startx = [fNewS0Start; fNewT2Start]; 
%     f = @(startx)OneCompSSD(startx, vY, vEchoTimes);
%     [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
%     vResnorms(i) = RESNORM; 
%     vParams(i, :) = vNonLinearParams; 
% end
% 
% %% Find the number of global minimums found
% vRoundedResnorm = round(vResnorms, 4); 
% nNumGlobalMins = length(find(vRoundedResnorm == min(vRoundedResnorm))); 
% disp(nNumGlobalMins);
% 
% 
% %% Using the fmincon to map the entire brain 
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonLinearS0Map = zeros(nDimX, nDimY, nDimZ);
mNonLinearT2Map = zeros(nDimX, nDimY, nDimZ); 
totalRESNORM = 0;
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                startx = [mNonNegS0Map(nRow, nCol, nSlice); mNonNegT2Map(nRow, nCol, nSlice)]; 
                [vParams, RESNORM] = NonLinearFit(vY, vEchoTime, startx, 3); 
                mNonLinearS0Map(nRow, nCol, nSlice) = vParams(1);
                mNonLinearT2Map(nRow, nCol, nSlice) = vParams(2);
                totalRESNORM = totalRESNORM + RESNORM;
            end
        end
    end
end
toc;


%% Functions 
function [vParams, fLowestResnorm] = NonLinearFit(vY, vEchoTimes, startx, iters)
    % 
    % Input 
    %      vY - n x 1 column vector of signal observations 
    %      vEchoTime - n x 1 column vector of corresponding echo time 
    %      startx - 2 x 1 column vector for the initial values of the non
    %               linear optmization 
    %      iters (int) - number of iterations per voxel to find the global
    %                    minimum
    % Output
    %      vParams - 2 x 1 column vector of the estimated parameters 
    % 
    % uses the fmincon to estimate parameters for a set of data
    % observations 
    
    A = [-1.0 0; 0 -1.0]; 
    b = [0; 0]; 
    f = @(startx)OneCompSSD(startx, vY, vEchoTimes);
    [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
    nIter = iters; 
    fLowestResnorm = RESNORM;
    vParams = vNonLinearParams; 
    for i = 1:nIter-1 
        fNewS0Start = startx(1) + randn()*500; 
        fNewT2Start = startx(2) + randn()*50;
        if(fNewS0Start < 0) 
            fNewS0Start = 1; 
        end
        if(fNewT2Start < 0)
            fNewT2Start = 1;
        end
        vNewStartx = [fNewS0Start; fNewT2Start]; 
        f = @(startx)OneCompSSD(startx, vY, vEchoTimes);
        [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, vNewStartx, A, b);
        if(RESNORM < fLowestResnorm)
            fLowestResnorm = RESNORM;
            vParams = vNonLinearParams;
        end
    end
end

function vParams = NonnegLeastSquaresFit(vY, vEchoTime)
    % Inputs:
    %      vY - n x 1 column vector of signal observations 
    %      vEchoTime - n x 1 column vector of corresponding echo time 
    % Output: 
    %      vParams - 2 x 1 column vector for the parameters S0 and T2
    % Uses the non negative Linear Least Squares method to estimate parameters 
    vNegativeIndices = find(vY <= 0); 
    while(length(vNegativeIndices) > 0)
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
    % corresponding echo time values.     vEchoTime = cell2mat(EchoTimes(nSubject)); 
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

function [fS0s, fT2s, bootStrapData] = BootStrap(vEchoTimes, S, startx, variance, h)
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
        [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fminunc('OneCompSSD', startx, h, bootStrapData(:, i), vEchoTimes); 
        disp(vNonLinearParams);
        fS0s(1, i) = vNonLinearParams(1);
        fT2s(1, i) = vNonLinearParams(2);
    end
end

function [fS0s, fT2s, bootStrapData] = BootStrap2(vEchoTimes, S, startx, variance)
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
        A = [-1 0; 0 -1];
        b = [0; 0];
        f = @(startx)OneCompSSD(startx, bootStrapData(:, i), vEchoTimes);
        [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b); 
        fS0s(1, i) = vNonLinearParams(1);
        fT2s(1, i) = vNonLinearParams(2);
    end
end