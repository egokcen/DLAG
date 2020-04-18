function gp_params = getGPparams_dlag(params, binWidth)
%
% gp_params = getGPparams_dlag(params, binWidth)
%
% Description: Delays and GP timescales can be found in params, but they
%              are easier to interpret when given in units of time. This
%              function gets within- and across-group timescales and delays
%              from params and converts them into the units of time
%              corresponding to binWidth. binWidth should match the 
%              binWidth to which the model was fitted.
%
% Arguments: 
%
%     params    -- Structure containing DLAG model parameters. 
%                  Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'stepSize ./ sqrt(gamma)'                                                    
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
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%
% Outputs:
%     
%    gp_params -- structure containing DLAG GP parameters, converted into
%                 units of time.
%                 DelayMatrix  -- (numGroups x xDim_across) array;
%                                 delays from across-group latents to 
%                                 observed variables
%                 tau_across -- (1 x xDim_across) array; across-group GP
%                               timescales
%                 tau_within -- (1 x numGroups) cell array; within-group
%                               GP timescales for each group. tau_within(i)
%                               is empty for groups with no within-group
%                               latents.
%                                     
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     10 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality

numGroups = length(params.xDim_within);

% Convert delays from bins to units of time
gp_params.DelayMatrix = binWidth .* params.DelayMatrix;

% Convert timescales to units of time
gp_params.tau_across = binWidth ./ sqrt(params.gamma_across);
gp_params.tau_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Leave tau_within empty for groups with no within-group latents
    if params.xDim_within(groupIdx) > 0
        gp_params.tau_within{groupIdx} = binWidth ./ sqrt(params.gamma_within{groupIdx});
    end
end