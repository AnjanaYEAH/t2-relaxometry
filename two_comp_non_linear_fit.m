%% Two compartment model fitting 

%% Two Compartment Non linear fitting of one compartment 
% We will use the non negative least squares as initial staring point for
% our iterative approach of finding the minimum 


%Define various options for the non-linear fitting algorithm 
h = optimset('MaxFunEvals', 2000, 'Algorithm', 'quasi-newton','Display', 'off', 'TolX', 1e-6,'TolFun', 1e-6);

%% Running the Two Comp Model once on a single voxel
%Get measurement mean and sigma 
[vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 57, 12, 29); 

%==========================================================================
% Starting parameter choices for Two Compartment 
%==========================================================================

vNonNegParams = NonnegLeastSquaresFit(vY, vEchoTimes); % Get starting point for non-linear fitting
startx = [sqrt(vNonNegParams(1)) 1/(0.5^2+1) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2))]; % Initial stating point 
[vBestParams, RESNORM] = TwoCompNonLinearFit(vY, vEchoTimes, startx, 10, h); % Two Compartment Non Linear Fitting 
S0 = vBestParams(1)^2;
v1 = 1/(vBestParams(2)^2+1);
T2_1 = vBestParams(3)^2;
T2_2 = vBestParams(4)^2;
v2 = 1-v1;
vSignal = S0*(v1*exp(-vEchoTime/T2_1) + v2*exp(-vEchoTime/T2_2));
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonLinearMylinMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearOccipetalGMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearOccipetalWMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearGMWMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearCSFMap = zeros(nDimX, nDimY, nDimZ);
totalRESNORM = 0;
mylincount = 0;
GMWMcount = 0;
CSFcount = 0;
mylinaccurate = 0;
CSFaccurate = 0;
GMWMaccurate = 0;
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vNonNegParams = NonnegLeastSquaresFit(vY, vEchoTimes); % Get starting point for non-linear fitting
                startx = [sqrt(vNonNegParams(1)) 1/(0.5^2+1) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2))]; % Initial stating point 
                [vBestParams, RESNORM] = TwoCompNonLinearFit(vY, vEchoTimes, startx, 3, h); % Two Compartment Non Linear Fitting
                S0 = vBestParams(1)^2;
                v1 = 1/(vBestParams(2)^2+1);
                T2_1 = vBestParams(3)^2;
                T2_2 = vBestParams(4)^2;
                v2 = 1-v1;
                totalRESNORM = totalRESNORM + RESNORM;
%                 %occipetal Grey matter
%                 if 41.6<T2_1 && T2_1<51.8
%                     mNonLinearOccipetalGMMap(nRow, nCol, nSlice) = T2_1;
%                 end
% %                 if 41.6<T2_2 && T2_2<51.8
% %                     mNonLinearOccipetalGMMap(nRow, nCol, nSlice) = T2_2;
% %                 end
%                 %occipetal WM matter
%                 if 44.7<T2_1 && T2_1<48.4
%                     mNonLinearOccipetalWMMap(nRow, nCol, nSlice) = T2_1;
%                 end
% %                 if 41.6<T2_2 && T2_2<51.8
% %                     mNonLinearOccipetalWMMap(nRow, nCol, nSlice) = T2_2;
% %                 end
%                 %Mylin
%                 if 1<T2_1 && T2_1<20
%                     mNonLinearMylinMap(nRow, nCol, nSlice) = T2_1;
%                     mylincount = mylincount + 1;
%                     if mCase01seg(nRow, nCol, nSlice, 4)
%                         mylinaccurate = mylinaccurate + 1;
%                     end
%                 end
%                 if 1<T2_2 && T2_2<20
%                     mNonLinearMylinMap(nRow, nCol, nSlice) = T2_2;
%                     mylincount = mylincount + 1;
%                     if mCase01seg(nRow, nCol, nSlice, 4)
%                         mylinaccurate = mylinaccurate + 1;
%                     end
%                 end
%                 %Grey matter and white matter
%                 if 80<T2_1 && T2_1<110
%                     mNonLinearGMWMMap(nRow, nCol, nSlice) = T2_1;
%                     GMWMcount = GMWMcount + 1;
%                     if mCase01seg(nRow, nCol, nSlice, 3) && mCase01seg(nRow, nCol, nSlice, 4)
%                         GMWMaccurate = GMWMaccurate + 1;
%                     end
%                 end
% %                 if 80<T2_2 && T2_2<110
% %                     mNonLinearGMWMMap(nRow, nCol, nSlice) = T2_2;
% %                 end
%                 %CSF 
%                 if 200<T2_1 && T2_1<600
%                     mNonLinearCSFMap(nRow, nCol, nSlice) = T2_1;
%                     CSFcount = CSFcount + 1;
%                     if mCase01seg(nRow, nCol, nSlice, 2)
%                         CSFaccurate = CSFaccurate + 1;
%                     end
%                 end
% %                 if 200<T2_2 && T2_2<600
% %                     mNonLinearCSFMap(nRow, nCol, nSlice) = T2_2;
% %                 end
            end
        end
    end
end
toc;
% imshow(flipud(squeeze(mNonLinearGMWMMap(:,:,30))'), []);
% [S0s, v1s, v2s, T2_1s, T2_2s, bootStrapData] = BootStrap(vEchoTimes, vSignal, vBestParams, h);
% histogram(S0s); 
% xlabel('S0');
% ylabel('Frequency');
% histogram(T2_1s); 
% xlabel('T2_1');
% ylabel('Frequency');
% histogram(T2_2s); 
% xlabel('T2_2');
% ylabel('Frequency');
% histogram(v1s); 
% xlabel('v1');
% ylabel('Frequency');
% histogram(v2s); 
% xlabel('v2');
% ylabel('Frequency');
% rangesS0s = ConfidenceInterval(S0s);
% rangesT2a= ConfidenceInterval(T2_1s);
% rangesT2b = ConfidenceInterval(T2_2s);
% rangesv1 = ConfidenceInterval(v1s);
% %% Getting the Two Compartment fitted parameters and synthetic singals 
% vBestParamsCorrected = [vBestParams(1)^2 1/(vBestParams(2)^2+1) vBestParams(3)^2 vBestParams(4)^2]; %Correct the contrained parameters
% vYEst = GetTwoCompSyntheticSignal(vBestParamsCorrected, vEchoTimes);
% fSSE = sum((vYEst - vY).^2, 'all');


%% Functions 

function [vParams, fLowestResnorm] = TwoCompNonLinearFit(vY, vEchoTimes, startx, iters, h)
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

    [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('TwoCompSSD',startx, h,vY, vEchoTimes);
    nIter = iters; 
    fLowestResnorm = RESNORM;
    vParams = vNonLinearParams; 
    for i = 1:nIter-1 
        fNewS0Start = normrnd(startx(1), abs(startx(1)));
        fNewVStart = normrnd(startx(2), abs(startx(2)));
        fNewT2aStart = normrnd(startx(3), abs(startx(3)));
        fNewT2bStart = normrnd(startx(4), abs(startx(4)));
        vNewStartx = [fNewS0Start; fNewVStart; fNewT2aStart; fNewT2bStart]; 
        [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('TwoCompSSD',vNewStartx,h,vY, vEchoTimes);
        if(RESNORM < fLowestResnorm)
            fLowestResnorm = RESNORM;
            vParams = vNonLinearParams;
        end
    end
end

function fLL = ComputeTwoCompLogLL(vY, dYSigma, vParams, vEchoTime)
    vSignals = GetTwoCompSyntheticSignal(vParams, vEchoTime);
    nN = length(vSignals); 
    dSSE = sum((vSignals - vY).^2);
    fLL = ((-nN/2) * log(2*pi*dYSigma^2)) - ((1/(2*dYSigma^2) * dSSE)); 
end


function vParams = NonLinearFit(vY, vEchoTimes, startx, iters)
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

function vEstY = GetTwoCompSyntheticSignal(vParams, vEchoTime)
    fS0 = vParams(1); 
    fVa = vParams(2); 
    fT2a = vParams(3);
    fT2b = vParams(4); 
    fVb = 1 - fVa;
    vEstY = fS0*(fVa*exp(-vEchoTime/fT2a) + fVb*exp(-vEchoTime/fT2b));
end 

function [S0s, v1s, v2s, T2_1s, T2_2s, bootStrapData] = BootStrap(vEchoTimes, S, startx, h)
    [nRows, nCols] = size(S);
    bootStrapSize = 500;
    bootStrapData = zeros(length(S), bootStrapSize);
    for i = 1: length(S)
        for j = 1: bootStrapSize
            bootStrapData(i, j) = S(i) + normrnd(0, sqrt(100));
        end
    end
    S0s = zeros(1, bootStrapSize);
    v1s = zeros(1, bootStrapSize);
    T2_1s = zeros(1, bootStrapSize);
    T2_2s = zeros(1, bootStrapSize);
    v2s = zeros(1, bootStrapSize);
    for i = 1: bootStrapSize
%         A = [-1.0 0 0 0; 0 -1.0 0 0; 0 0 -1.0 0; 0 0 0 -1.0]; 
%         b = [0; 0; 0; 0]; 
%         f = @(startx)TwoCompSSD(startx, bootStrapData(:, i), vEchoTimes);
%         [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
        [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('TwoCompSSD',startx,h,bootStrapData(:, i), vEchoTimes);
        vBestParams = vNonLinearParams;
        S0s(1, i) = vBestParams(1)^2;
        v1s(1, i) = 1/(vBestParams(2)^2+1);
        T2_1s(1, i) = vBestParams(3)^2;
        T2_2s(1, i) = vBestParams(4)^2;
        v2s(1, i) = 1-v1s(1, i);
    end
end