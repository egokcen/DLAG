function [varexp, domexp] = computeVarExp_dlag(params, total)
%
% [varexp, domexp] = computeVarExp_dlag(params)
%
% Description: Compute the fraction of variance explained in each group by 
%              latent variables. Choose the denominator of the computation
%              (total variance or shared variance) via the 'total'
%              argument.
%
% Arguments:
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    d            -- (yDim x 1) array; observation mean
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    R            -- (yDim x yDim) array; observation noise
%                                    covariance 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     total -- logical; set true to set total variance as the denominator,
%              false to set shared variance as the denominator.
%
% Outputs:
%
%     varexp -- structure with the following fields:
%                   across -- (1 x numGroups) array; fraction variance 
%                             explained by across-group latents.
%                   within -- (1 x numGroups) array; fraction variance 
%                             explained by within-group latents.
%                   indiv  -- (1 x numGroups) cell array; indiv{i} is an
%                             array with the variance explained by each
%                             individual latent variable.
%                   dom    -- (1 x numGroups) cell array; dom{i} gives the
%                             variance explained by each dominant mode.
%     domexp -- structure with the following fields:
%                   across -- (1 x numGroups) array; 
%                             across{i} -- (xDim(i) x 1) array; fraction 
%                             of each dominant mode explained by
%                             across-group latents.
%                   within -- (1 x numGroups) array; 
%                             within{i} -- (xDim(i) x 1) array; fraction 
%                             of each dominant mode explained by
%                             within-group latents.
%                   indiv  -- (1 x numGroups) cell array; 
%                             indiv{i} -- (xDim(i) x xDim(i)); each row 
%                             gives the fraction of a dominant mode 
%                             explained by each individual latent variable.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Mar 2021 -- Initial full revision.
%     18 Mar 2021 -- Added fraction of dominant modes explained.

% Initialize output structure
numGroups = length(params.yDims);
varexp.across = nan(1,numGroups);
varexp.within = nan(1,numGroups);
varexp.indiv = cell(1,numGroups);
varexp.dom = cell(1,numGroups);
domexp.across = cell(1,numGroups);
domexp.within = cell(1,numGroups);
domexp.indiv = cell(1,numGroups);

% Compute dominant modes, including within- and across-group latents
[S, ~, ~, H] = dominantModes_dlag(params);

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);

% Compute variance explained
for groupIdx = 1:numGroups
    
    % All parameters for the current group
    R = groupParams{groupIdx}.R;
    C = groupParams{groupIdx}.C;
    Cvars = vecnorm(C,2,1).^2; % Variances of individual dimensions
    
    % Across-group parameters
    acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:params.xDim_across, cell(1,numGroups));
    Ca = acrossParams.C;
    
    % Within-group parameters
    withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:params.xDim_within(groupIdx)});
    Cw = withinParams.C;
    
    % Compute fraction of variance explained
    denom = 0; % Denominator of computation
    if total
        % Use total variance 
        denom = trace(Ca*Ca' + Cw*Cw' + R);
    else
        % Use shared variance
        denom = trace(Ca*Ca' + Cw*Cw');
    end
    varexp.across(groupIdx) = trace(Ca*Ca') / denom;
    varexp.within(groupIdx) = trace(Cw*Cw') / denom;
    varexp.indiv{groupIdx} = Cvars ./ denom;
    varexp.dom{groupIdx} = (diag(S{groupIdx}).^2)' ./ denom;
    
    % Fraction of dominant modes explained
    domexp.indiv{groupIdx} = (H{groupIdx}.^2) ./ (diag(S{groupIdx}).^2);
    domexp.across{groupIdx} = sum(domexp.indiv{groupIdx}(:,1:params.xDim_across),2);
    domexp.within{groupIdx} = sum(domexp.indiv{groupIdx}(:,params.xDim_across+1:end),2);
    
end