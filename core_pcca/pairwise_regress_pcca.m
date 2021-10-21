function [Ypred, R2, MSE] = pairwise_regress_pcca(Ys, params, rGroups)
%
% [Ypred, R2, MSE] = pairwise_regress_pcca(Ys, params, rGroups)
%
% Description: Performs regression between two groups using a pCCA model.
%
% Arguments:
%     Ys     -- (1 x numGroups) cell array; list of data matrices 
%               {(y1Dim x N), (y2Dim x N), ...}
%     params -- learned PCCA parameteres (structure with fields Cs, Rs, ds)
%     rGroups -- (1 x 2) array; Each element specifies a group to be 
%                included in the regression.
%
% Outputs:
%     Ypred  -- (1 x 2) cell array; predicted responses of target group
%               given responses of source group
%     R2     -- (1 x 2) array; R^2 in each pairwise direction
%     MSE    -- (1 x 2) array; mean-squared error in each pairwise direction
%
%     NOTE: The group order of outputs corresponds to the group order of
%           rGroups. For example, R2(1) is the performance of group 
%           rGroups(1) predicting group rGroups(2). R2(2) is the 
%           performance of group rGroups(2) predicting group rGroups(1).
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.
%     20 Oct 2021 -- Updated documentation to clarify group order of
%                    outputs.

numRGroups = length(rGroups); % Really just 2, for pairwise regression
Ypred = cell(1,numRGroups);
R2 = nan(1,numRGroups);
MSE = nan(1,numRGroups);
for rIdx = 1:numRGroups
    in = rGroups(rIdx); % Source group
    out = rGroups(setdiff(1:numRGroups,rIdx)); % Target group
    
    Yin = Ys{in};
    Cin = params.Cs{in};
    Rin = params.Rs{in};
    din = params.ds{in};
    
    Yout = Ys{out};
    Cout = params.Cs{out};
    dout = params.ds{out};
    
    N = size(Yin,2);
    
    % Predict
    Ypred{rIdx} = Cout * Cin' / (Cin * Cin' + Rin) * (Yin - repmat(din,1,N)) ...
                + repmat(dout,1,N);
    
    % Evaluate performance
    RSS = sum( sum( ( Yout - Ypred{rIdx} ).^2, 1 ) );
    TSS = sum( sum( ( Yout - repmat( mean(Yout,2), [1 size(Yout,2)] ) ).^2, 1 ) );
    R2(rIdx) = 1 - RSS / TSS;
    MSE(rIdx) = immse(Ypred{rIdx}, Yout);

end