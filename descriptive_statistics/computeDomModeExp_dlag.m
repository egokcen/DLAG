function domexp = computeDomModeExp_dlag(params)
%
% domexp = computeDomModeExp_dlag(params, total)
%
% Description: Compute the fraction of each dominant mode explained by 
%              latent variables.
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
% Outputs:
%
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
%     18 Mar 2021 -- Initial full revision.

% Initialize output structure
numGroups = length(params.yDims);
domexp.across = cell(1,numGroups);
domexp.within = cell(1,numGroups);
domexp.indiv = cell(1,numGroups);

% Compute dominant modes, including within- and across-group latents
[S, ~, ~, H] = dominantModes_dlag(params);

% Compute fraction of each dominant mode explained by latents
for groupIdx = 1:numGroups
    
    
    domexp.across(groupIdx) = trace(Ca*Ca') / denom;
    domexp.within(groupIdx) = trace(Cw*Cw') / denom;
    domexp.indiv{groupIdx} = Cvars ./ denom;
end