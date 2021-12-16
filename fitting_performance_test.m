%% This File produces data on the time cost of the fitting algorithms
% We save the time elapsed for fitting a single voxel. We then fit across
% slice 30 to get a sequence of time values. 

%% Linear Fitting Time Cost 
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mLinearS0Map = zeros(nDimX, nDimY, nDimZ);
mLinearT2Map = zeros(nDimX, nDimY, nDimZ);
% Stores the time computation of linear least squares 
vLinearTimeElapsed = zeros(nDimX * nDimY, 1); 
nLinearTimeElapsedIndex = 0; 
% Stores the residual of the linear least squares 
vLinearRes = zeros(nDimX * nDimY, 1); 
vLinearResC = zeros(nDimX * nDimY, 1);
nLinearResIndex = 0; 
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = LinearLeastSquaresFit(vY, vEchoTime); 
                %tempText = sprintf("Row: %d | Col: %d", nRow, nCol);
                %disp(tempText);
                %disp(vParams);
                mLinearS0Map(nRow, nCol, nSlice) = vParams(1);
                mLinearT2Map(nRow, nCol, nSlice) = vParams(2);
                % Store the time elapsed 
                nLinearTimeElapsedIndex = nLinearTimeElapsedIndex + 1; 
                vLinearTimeElapsed(nLinearTimeElapsedIndex) = toc; 
                % Store the residual 
                nLinearResIndex = nLinearResIndex + 1; 
                vLinearRes(nLinearResIndex)= OneCompSSD([vParams(1) vParams(2)], vY, vEchoTime);
                if(nLinearResIndex == 1) 
                    vLinearResC(1) = vLinearRes(1);
                else
                    vLinearResC(nLinearResIndex) = vLinearRes(nLinearResIndex) + vLinearResC(nLinearResIndex-1);
                end
            end
        end
    end
end

%% Display Linear Fitting time plot
figure(1); 
subplot(1, 2, 1); 
plot(linspace(1, nLinearTimeElapsedIndex, nLinearTimeElapsedIndex), vLinearTimeElapsed(1:nLinearTimeElapsedIndex));
title("Linear Fit Time Elapsed against Voxels"); 
xlabel("Number of voxels"); 
ylabel("Seconds");

%% Display Linear Fitting Res plot per 
subplot(1, 2, 2); 
%bar(linspace(1, nLinearResIndex, nLinearResIndex), vLinearRes(1:nLinearResIndex)); 
plot(linspace(1, nLinearResIndex, nLinearResIndex), vLinearResC(1:nLinearResIndex));
title("Linear Fit Cumulative Residual"); 
xlabel("Voxel number"); 
ylabel("Residual"); 

%% Non Negative Linear Fitting Time Cost
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonNegLinearS0Map = zeros(nDimX, nDimY, nDimZ);
mNonNegLinearT2Map = zeros(nDimX, nDimY, nDimZ);
% Stores the time computation for the Non Negative Linear fitting 
vNonNegLinearTimeElapsed = zeros(nDimX * nDimY, 1); 
nNonNegLinearTimeElapsedIndex = 0; 
% Store the residual for the Non Negative Linear fitting 
vNonNegLinearRes = zeros(nDimX * nDimY, 1); 
vNonNegLinearResC = zeros(nDimX * nDimY, 1); 
nNonNegLinearResIndex = 0; 
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = NonnegLeastSquaresFit(vY, vEchoTime); 
                mNonNegLinearS0Map(nRow, nCol, nSlice) = vParams(1);
                mNonNegLinearT2Map(nRow, nCol, nSlice) = vParams(2);
                % Store the time elapsed 
                nNonNegLinearTimeElapsedIndex = nNonNegLinearTimeElapsedIndex + 1; 
                vNonNegLinearTimeElapsed(nNonNegLinearTimeElapsedIndex) = toc; 
                % Store the residual of this parameter estimate 
                nNonNegLinearResIndex = nNonNegLinearResIndex + 1; 
                vNonNegLinearRes(nNonNegLinearResIndex)= OneCompSSD([vParams(1) vParams(2)], vY, vEchoTime);
                if(nNonNegLinearResIndex == 1) 
                    vNonNegLinearResC(1) = vNonNegLinearRes(1);
                else
                    vNonNegLinearResC(nNonNegLinearResIndex) = vNonNegLinearRes(nNonNegLinearResIndex) + vNonNegLinearResC(nNonNegLinearResIndex-1);
                end
            end
        end
    end
end


%% Display Non Negative Linear Fitting time plot
figure(2); 
subplot(1, 2, 1);
plot(linspace(1, nNonNegLinearTimeElapsedIndex, nNonNegLinearTimeElapsedIndex), vNonNegLinearTimeElapsed(1:nNonNegLinearTimeElapsedIndex));
title("Non Negative Linear Fit Time Elapsed against Voxels"); 
xlabel("Number of voxels"); 
ylabel("Seconds");

%% Display Non Negative Linear Fitting Res plot per 
subplot(1, 2, 2); 
plot(linspace(1, nNonNegLinearResIndex, nNonNegLinearResIndex), vNonNegLinearResC(1:nNonNegLinearResIndex));
title("NonNeg Linear Fit Cumulative Residual"); 
xlabel("Voxel number"); 
ylabel("Residual"); 

%% Weighted least squares Fitting time cost
%
% In this section we encapsulate finding the parameters for a single voxel
% into a single function called LinearLeastSquares() takes ~ 90 seconds
nDimX = 96;
nDimY = 96;
nDimZ = 55;
mWS0Map = zeros(nDimX, nDimY, nDimZ);
mWT2Map = zeros(nDimX, nDimY, nDimZ);
vWeightedTimeElapsed = zeros(nDimX * nDimY, 1); 
nWeightedTimeElapsedIndex = 0; 
% Stores the residual of the linear least squares 
vWLinearRes = zeros(nDimX * nDimY, 1); 
vWLinearResC = zeros(nDimX * nDimY, 1);
nWLinearResIndex = 0; 
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                %This voxel will be mapped 
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                vParams = WeightedLeastSquaresFit(vY, vEchoTime); 
                mWS0Map(nRow, nCol, nSlice) = vParams(1);
                mWT2Map(nRow, nCol, nSlice) = vParams(2);
                nWeightedTimeElapsedIndex = nWeightedTimeElapsedIndex + 1; 
                vWeightedTimeElapsed(nWeightedTimeElapsedIndex) = toc;
                % Store the residual 
                nWLinearResIndex = nWLinearResIndex + 1; 
                vWLinearRes(nWLinearResIndex)= OneCompSSD(vParams, vY, vEchoTime);
                if(nWLinearResIndex == 1) 
                    vWLinearResC(1) = vWLinearRes(1);
                else
                    vWLinearResC(nWLinearResIndex) = vWLinearRes(nWLinearResIndex) + vWLinearResC(nWLinearResIndex-1);
                end
            end
        end
    end
end

%% Display Weighted Linear Fitting time plot
figure(3); 
subplot(1, 2, 1);
plot(linspace(1, nWeightedTimeElapsedIndex, nWeightedTimeElapsedIndex), vWeightedTimeElapsed(1:nWeightedTimeElapsedIndex));
title("Weighted Linear Fit Time Elapsed against Voxels"); 
xlabel("Number of voxels"); 
ylabel("Seconds");

%% Display Weighted Lienar Fitting Res plot per 
subplot(1, 2, 2); 
plot(linspace(1, nWLinearResIndex, nWLinearResIndex), vWLinearResC(1:nWLinearResIndex));
title("Weighted Linear Fit Cumulative Residual"); 
xlabel("Voxel number"); 
ylabel("Residual"); 


%% Display Linear, Non Negative Linear and Weighted time graphs into single figure 
figure(4); 
subplot(1, 2, 1);
plot(linspace(1, nLinearTimeElapsedIndex, nLinearTimeElapsedIndex), vLinearTimeElapsed(1:nLinearTimeElapsedIndex));
hold on 
plot(linspace(1, nNonNegLinearTimeElapsedIndex, nNonNegLinearTimeElapsedIndex), vNonNegLinearTimeElapsed(1:nNonNegLinearTimeElapsedIndex));
plot(linspace(1, nWeightedTimeElapsedIndex, nWeightedTimeElapsedIndex), vWeightedTimeElapsed(1:nWeightedTimeElapsedIndex));
hold off
title("Linear fit time computation cost");
xlabel("Number of voxels"); 
ylabel("Seconds");
legend('Linear LeastSqr', 'NonNeg Linear LeastSqr', 'Weighted Least Squares');

nUpTo = 150;
subplot(1, 2, 2); 
p = plot(linspace(1, nUpTo, nUpTo), vLinearResC(1:nUpTo));
p(1).LineStyle = '-';
p(1).LineWidth = 2;
hold on
plot(linspace(1, nUpTo, nUpTo), vNonNegLinearResC(1:nUpTo), 'LineWidth', 1);
plot(linspace(1, nUpTo, nUpTo), vWLinearResC(1:nUpTo));
hold off
title("Linear fit Cumulative Residuals"); 
xlabel("Voxel number"); 
ylabel("Residual"); 
legend('Linear LeastSqr', 'NonNeg Linear LeastSqr', 'Weighted Least Squares');
%% Non Linear Fitting Time Cost 
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonLinearS0Map = zeros(nDimX, nDimY, nDimZ);
mNonLinearT2Map = zeros(nDimX, nDimY, nDimZ); 
vNonLinearTimeElapsed = zeros(nDimX * nDimY, 1); 
nNonLinearTimeElapsedIndex = 0; 
tic;
for nRow = 1:nDimX
    disp(nRow);
    for nCol = 1:nDimX
        for nSlice = 30:30
            if(GetGlobalMaskVal(Masks, nRow, nCol, nSlice))
                [vY, vEchoTime] = GetVoxelMeasurements(Data, ET, nRow, nCol, nSlice);
                startx = [mNonNegS0Map(nRow, nCol, nSlice); mNonNegT2Map(nRow, nCol, nSlice)]; 
                vParams = NonLinearFit(vY, vEchoTime, startx, 1);  %The last parameter here set how many iterations we do per search
                mNonLinearS0Map(nRow, nCol, nSlice) = vParams(1);
                mNonLinearT2Map(nRow, nCol, nSlice) = vParams(2);
                nNonLinearTimeElapsedIndex = nNonLinearTimeElapsedIndex + 1; 
                vNonLinearTimeElapsed(nNonLinearTimeElapsedIndex) = toc; 
            end
        end
    end
end

%% Display Non Linear Fitting time plot
figure(5); 
plot(linspace(1, nNonLinearTimeElapsedIndex, nNonLinearTimeElapsedIndex), vNonLinearTimeElapsed(1:nNonLinearTimeElapsedIndex));
title("Non Linear Fitting time cost per voxel");
xlabel("Number of voxels"); 
ylabel("Seconds");

%% Display all Fitting time plots together 
figure(6); 
nUpTo = 5;
plot(linspace(1, nUpTo, nUpTo), vLinearTimeElapsed(1:nUpTo));
hold on 
plot(linspace(1, nUpTo, nUpTo), vNonNegLinearTimeElapsed(1:nUpTo));
plot(linspace(1, nUpTo, nUpTo), vNonLinearTimeElapsed(1:nUpTo));
hold off
xlabel("Number of voxels"); 
ylabel("Seconds");
legend('Linear LeastSqr', 'NonNeg Linear LeastSqr', 'Non-linear LeastSqr');


%% Functions 
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

