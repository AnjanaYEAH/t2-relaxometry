h = optimset('MaxFunEvals', 2000, 'Algorithm', 'quasi-newton','Display', 'off', 'TolX', 1e-6,'TolFun', 1e-6);
[vY, vEchoTimes] = GetVoxelMeasurements(Data, ET, 57, 12, 29); 
r = GetHighestNfractions(Dataseg, 57, 12, 29, 3);
vNonNegParams = NonnegLeastSquaresFit(vY, vEchoTimes); % Get starting point for non-linear fitting
startx = [sqrt(vNonNegParams(1)) asin(sqrt(r(1))) asin(sqrt(r(2))) asin(sqrt(r(3))) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2))]; % Initial stating point
% startx = [sqrt(vNonNegParams(1)) 1/(sqrt(0.5)-1) 1/(sqrt(0.5)-2) 1/(sqrt(0.5)-3) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2))]; % Initial stating point 

[vBestParams, RESNORM] = ThreeCompNonLinearFit(vY, vEchoTimes, startx, 10, h); % Two Compartment Non Linear Fitting 
S0 = vBestParams(1)^2;
v1 = sin(vBestParams(2))^2;
v2 = sin(vBestParams(3))^2;
v3 = sin(vBestParams(4))^2;
a = [v1, v2, v3];
a = sort(a);
total = v1 + v2 + v3;
v3 = a(1)/total;
v2 = a(2)/total;
v1 = a(3)/total;
T2_1 = vBestParams(5)^2;
T2_2 = vBestParams(6)^2;
T2_3 = vBestParams(7)^2;
vSignal = S0*(v1*exp(-vEchoTime/T2_1) + v2*exp(-vEchoTime/T2_2) + v3*exp(-vEchoTime/T2_3));
[S0s, v1s, v2s, v3s, T2_1s, T2_2s, T2_3s, bootStrapData] = BootStrap(vEchoTimes, vSignal, vBestParams, h);

% histogram(S0s);
% xlabel('S0');
% ylabel('Frequency');
% rangesS0s = ConfidenceInterval(S0s);
% 
% histogram(v1s);
% xlabel('v1');
% ylabel('Frequency');
% rangesv1s = ConfidenceInterval(v1s);
% 
% histogram(v2s);
% xlabel('v2');
% ylabel('Frequency');
% rangesv2s = ConfidenceInterval(v2s);
%  
% histogram(v3s);
% xlabel('v3');
% ylabel('Frequency');
% rangesv3s = ConfidenceInterval(v3s);
%  
% histogram(T2_1s);
% xlabel('T2_1');
% ylabel('Frequency');
% rangesT2_1 = ConfidenceInterval(T2_1s);
% 
% histogram(T2_2s);
% xlabel('T2_2');
% ylabel('Frequency');
% rangesT2_2 = ConfidenceInterval(T2_2s);
% 
% histogram(T2_3s);
% xlabel('T2_3');
% ylabel('Frequency');
% rangesT2_3 = ConfidenceInterval(T2_3s);
% 
nDimX = 96;
nDimY = 96;
nDimZ = 55; 
mNonLinearMylinMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearOccipetalGMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearOccipetalWMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearGMWMMap = zeros(nDimX, nDimY, nDimZ);
mNonLinearCSFMap = zeros(nDimX, nDimY, nDimZ);
mylincount = 0;
GMWMcount = 0;
CSFcount = 0;
totalRESNORM = 0;
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
                r = GetHighestNfractions(Dataseg, nRow, nCol, nSlice, 3);
                vNonNegParams = NonnegLeastSquaresFit(vY, vEchoTimes); % Get starting point for non-linear fitting
                startx = [sqrt(vNonNegParams(1)) asin(sqrt(r(1))) asin(sqrt(r(2))) asin(sqrt(r(3))) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2)) sqrt(vNonNegParams(2))];
                [vBestParams, RESNORM] = ThreeCompNonLinearFit(vY, vEchoTimes, startx, 3, h); % Two Compartment Non Linear Fitting 
                S0 = vBestParams(1)^2;
                v1 = sin(vBestParams(2))^2;
                v2 = sin(vBestParams(3))^2;
                v3 = sin(vBestParams(4))^2;
                a = [v1, v2, v3];
                a = sort(a);
                total = v1 + v2 + v3;
                v3 = a(1)/total;
                v2 = a(2)/total;
                v1 = a(3)/total;
                T2_1 = vBestParams(5)^2;
                T2_2 = vBestParams(6)^2;
                T2_3 = vBestParams(7)^2;
                totalRESNORM = totalRESNORM + RESNORM;
%                 %occipetal Grey matter
%                 if 41.6<T2_1 && T2_1<51.8
%                     mNonLinearOccipetalGMMap(nRow, nCol, nSlice) = T2_1;
%                 end
% %                 if 41.6<T2_2 && T2_2<51.8
% %                     mNonLinearOccipetalGMMap(nRow, nCol, nSlice) = T2_2;
% %                 end
% %                 if 41.6<T2_3 && T2_3<51.8
% %                     mNonLinearOccipetalGMMap(nRow, nCol, nSlice) = T2_3;
% %                 end
%                 %occipetal WM matter
%                 if 44.7<T2_1 && T2_1<48.4
%                     mNonLinearOccipetalWMMap(nRow, nCol, nSlice) = T2_1;
%                 end
% %                 if 41.6<T2_2 && T2_2<51.8
% %                     mNonLinearOccipetalWMMap(nRow, nCol, nSlice) = T2_2;
% %                 end
% %                 if 41.6<T2_3 && T2_3<51.8
% %                     mNonLinearOccipetalWMMap(nRow, nCol, nSlice) = T2_3;
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
%                 if 1<T2_3 && T2_3<20
%                     mNonLinearMylinMap(nRow, nCol, nSlice) = T2_3;
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
% %                 if 80<T2_3 && T2_3<110
% %                     mNonLinearGMWMMap(nRow, nCol, nSlice) = T2_3;
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
% %                 if 200<T2_3 && T2_3<600
% %                     mNonLinearCSFMap(nRow, nCol, nSlice) = T2_3;
% %                 end
            end
        end
    end
end
toc;
% imshow(flipud(squeeze(mNonLinearMylinMap(:,:,30))'), []);
% countCSF = 0;
% countMyelin = 0;
% countGMWM = 0;
% for i = 1: 96
%     for j = 1: 96
%         if mCase01seg(i, j, 29, 2)
%             countCSF = countCSF + 1;
%         end
%         if mCase01seg(i, j, 29, 3) && mCase01seg(i, j, 29, 4)
%             countGMWM = countGMWM + 1;
%         end
%         if mCase01seg(i, j, 29, 4)
%             countMyelin = countMyelin + 1;
%         end
%     end
% end
function [vParams, fLowestResnorm] = ThreeCompNonLinearFit(vY, vEchoTimes, startx, iters, h)
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
    
%     A = [-1.0 0 0 0 0 0 0; 0 -1.0 0 0 0 0 0; 0 0 -1.0 0 0 0 0; 0 0 0 -1.0 0 0 0; 0 0 0 0 -1.0 0 0; 0 0 0 0 0 -1.0 0; 0 0 0 0 0 0 -1.0]; 
%     b = [0; 0; 0; 0; 0; 0; 0]; 
%     f = @(startx)ThreeCompSSD(startx, vY, vEchoTimes);
%     [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
    [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('ThreeCompSSD',startx,h,vY, vEchoTimes);
    nIter = iters; 
    fLowestResnorm = RESNORM;
    vParams = vNonLinearParams; 
    for i = 1:nIter-1 
        fNewS0Start = normrnd(startx(1), abs(startx(1)));
        fNewV1Start = normrnd(startx(2), abs(startx(2)));
        fNewV2Start = normrnd(startx(3), abs(startx(3)));
        fNewV3Start = normrnd(startx(4), abs(startx(4)));
        fNewT2aStart = normrnd(startx(5), abs(startx(5)));
        fNewT2bStart = normrnd(startx(6), abs(startx(6))); 
        fNewT2cStart = normrnd(startx(7), abs(startx(7)));
        vNewStartx = [fNewS0Start; fNewV1Start; fNewV2Start; fNewV3Start; fNewT2aStart; fNewT2bStart; fNewT2cStart]; 
        [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('ThreeCompSSD',vNewStartx,h,vY, vEchoTimes);
        if(RESNORM < fLowestResnorm)
            fLowestResnorm = RESNORM;
            vParams = vNonLinearParams;
        end
    end
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

function r = GetHighestNfractions(Dataseg, nX, nY, nZ, n)
    [dim1, dim2] = size(Dataseg);
    segdata = zeros(dim2, dim2);
    seg = zeros(dim1, dim2);
    for i =1:dim2
        d = Dataseg(dim1, i);
        d = cell2mat(d);
        segdata(i, :) = d(nX, nY, nZ, :);
    end 
    for i =1:dim2
        seg(1, i) = mean(segdata(:, i));
    end
    r = maxk(seg, n);
    for i =1:n
        if r(i) == 0
            r(i) = 0.1;
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

function [S0s, v1s, v2s, v3s, T2_1s, T2_2s, T2_3s, bootStrapData] = BootStrap(vEchoTimes, S, startx, h)
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
    v2s = zeros(1, bootStrapSize);
    v3s = zeros(1, bootStrapSize);
    T2_1s = zeros(1, bootStrapSize);
    T2_2s = zeros(1, bootStrapSize);
    T2_3s = zeros(1, bootStrapSize);
    for i = 1: bootStrapSize
%         A = [-1.0 0 0 0; 0 -1.0 0 0; 0 0 -1.0 0; 0 0 0 -1.0]; 
%         b = [0; 0; 0; 0]; 
%         f = @(startx)TwoCompSSD(startx, bootStrapData(:, i), vEchoTimes);
        [vNonLinearParams,RESNORM,EXITFLAG,OUTPUT]=fminunc('ThreeCompSSD',startx,h,bootStrapData(:, i), vEchoTimes);
%         [vNonLinearParams, RESNORM, EXITFLAG, OUTPUT] = fmincon(f, startx, A, b);
        vBestParams = vNonLinearParams;
        S0s(1, i) = vBestParams(1)^2;
        v1 = sin(vBestParams(2))^2;
        v2 = sin(vBestParams(3))^2;
        v3 = sin(vBestParams(4))^2;
        a = [v1, v2, v3];
        a = sort(a);
        total = v1 + v2 + v3;
        v1s(1, i) = a(3)/total;
        v2s(1, i) = a(2)/total;
        v3s(1, i) = a(1)/total;
        T2_1s(1, i) = vBestParams(5)^2;
        T2_2s(1, i) = vBestParams(6)^2;
        T2_3s(1, i) = vBestParams(7)^2;
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