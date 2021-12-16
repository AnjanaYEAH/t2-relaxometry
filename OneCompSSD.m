function fSumRes = OneCompSSD(vParams, vY, vEchoTime)
%ONECOMPSSD is the objective function for the one compartment model 
%   Detailed explanation goes here
    fS0 = vParams(1); 
    fT2 = vParams(2); 
    [nRows, nCols] = size(vY); 
    vEstY = zeros(nRows, 1); 
    
    for nRow = 1:nRows
        vEstY(nRow) = fS0 * exp(vEchoTime(nRow) * (-1.0 / fT2));
    end
    fSumRes = sum((vEstY - vY).^2); 
end

