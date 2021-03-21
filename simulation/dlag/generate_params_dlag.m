function params = generate_params_dlag(yDims, xDim_across, xDim_within, ...
    binWidth, snr, tau_lim, eps_lim, delay_lim)
% params = generate_params_dlag(...)
%
% Description: Randomly generate DLAG model parameters, within 
%              specified constraints. The params output structure is
%              compatible with other DLAG codepack functions.
%
% Arguments:
%
%     Required:
%
%     yDims       -- (1 x numGroups) array; List of dimensionalities of
%                    observed data, [Y1, Y2,...]
%     xDim_across -- int; Number of across-group latents
%     xDim_within -- (1 x numGroups) array; Number of within-group latents
%                    for each group. 
%     binWidth    -- float; intended spike count bin width or sample period 
%                    (in units of time). Assume uniform sampling.
%     snr         -- (1 x numGroups) array; List of signal-to-noise ratios,
%                    defined as trace(CC') / trace(R)
%     tau_lim   -- (1 x 2) array; lower- and upper-bounds of GP timescales
%     eps_lim   -- (1 x 2) array; lower- and upper-bounds of GP noise variances
%     delay_lim -- (1 x 2) array; lower- and upper-bounds of delays, in
%                  units of time
%
% Outputs:
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'binWidth ./ sqrt(gamma)'                                                    
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
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     08 Jan 2021 -- Initial full revision.

    % Constants
    numGroups = length(yDims);
    centered = false; % Allow for a non-zero mean parameter
    
    % Start filling in some basic fields in the output structure
    params.covType = 'rbf';
    params.xDim_across = xDim_across;
    params.xDim_within = xDim_within;
    params.yDims = yDims;
    
    % Generate GP parameters
    min_tau = tau_lim(1);
    max_tau = tau_lim(2);
    min_eps = eps_lim(1);
    max_eps = eps_lim(2);
    min_delay = delay_lim(1);
    max_delay = delay_lim(2);

    % Across-group timescales and noise variances
    taus_across = [];
    eps_across = [];
    if xDim_across > 0
        % If xDim_across is 0, keep taus and eps empty
        taus_across = min_tau + (max_tau-min_tau).*rand(1,xDim_across); 
        % Deal with noise variances on a log scale
        eps_across = exp(log(min_eps) + (log(max_eps)- log(min_eps)).*(rand(1,xDim_across)));
    end
    % Fill in output structure
    params.gamma_across = (binWidth ./ taus_across).^2;
    params.eps_across = eps_across;

    % Within-group timescales and noise variances
    taus_within = cell(1,numGroups);
    eps_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        % If xDim_within(groupIdx) is 0, keep taus and eps empty
        if xDim_within(groupIdx) > 0
            taus_within{groupIdx} = min_tau + (max_tau-min_tau).*rand(1,xDim_within(groupIdx));
            eps_within{groupIdx} = exp(log(min_eps) + (log(max_eps)- log(min_eps)).*(rand(1,xDim_within(groupIdx))));
        end
    end
    % Fill in output structure
    for groupIdx = 1:numGroups
        params.gamma_within{groupIdx} = (binWidth ./ taus_within{groupIdx}).^2;
        params.eps_within{groupIdx} = eps_within{groupIdx};
    end

    % Delays
    delays = cell(1,numGroups);
    if xDim_across > 0
        for groupIdx = 1:numGroups
            if groupIdx <= 1
                % Set delays to first group to 0 time steps
                delays{groupIdx} = zeros(1,xDim_across); 
            else
                % All other groups have non-zero delays
                delays{groupIdx} = min_delay + (max_delay-min_delay).*rand(1,xDim_across);
            end
        end
    end
    % Fill in output structure
    params.DelayMatrix = cat(1, delays{:})./binWidth;
    
    % Generate observation model parameters
    Cs = cell(1,numGroups);
    Rs = cell(1,numGroups);
    ds = cell(1,numGroups);
    for groupIdx = 1:numGroups
        % Total latent dimensionality of current group
        xDim = xDim_across + xDim_within(groupIdx);
        [C, R, d] = generate_pcca_params(yDims(groupIdx), xDim, snr(groupIdx), centered);
        % Collect observation parameters for the current group
        Cs{groupIdx} = C{1};
        % Take only the diagonal elements of R (pCCA generates a full
        % covariance matrix for each group)
        Rs{groupIdx} = diag(diag(R{1})); 
        ds{groupIdx} = d{1};
    end
    
    % Fill in output structure
    params.C = blkdiag(Cs{:});
    params.R = blkdiag(Rs{:});
    params.d = cat(1, ds{:});
    
end