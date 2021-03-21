function outparams = convertParamUnits(inparams, inres, outres)
%
% outparams = convertParamUnits(inparams, inres, outres)
%
% Description: Convert the (time) units of DLAG parameters based on sample
%              period 'inres' to units based on sample period 'outres'.
%
% Arguments:
%
%     inparams  -- Structure containing DLAG model parameters, in units
%                  normalized by input resolution, inres. 
%                  Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'inres ./ sqrt(gamma)'                                                    
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
%                                    time-steps, spaced according to inres.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     inres     -- float; input resolution (sample period or bin width), in
%                  units of time.
%     outres    -- float; output resolution, in units of time.
%
% Outputs:
%
%     outparams  -- Structure containing DLAG model parameters, in units
%                   normalized by output resolution, outres. 
%                   Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'outres ./ sqrt(gamma)'                                                    
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
%                                    time-steps, spaced according to outres.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Mar 2021 -- Initial full revision.  

numGroups = length(inparams.yDims);
outparams = inparams;
outparams.gamma_across = (outres / inres)^2 .* inparams.gamma_across;
outparams.DelayMatrix = (inres / outres) .* inparams.DelayMatrix;
for groupIdx = 1:numGroups
    outparams.gamma_within{groupIdx} = (outres / inres)^2 .* inparams.gamma_within{groupIdx}; 
end