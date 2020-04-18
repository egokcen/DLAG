function [Ys, Xs, params] = simdata_dlag(N, T, samplePeriod, yDims, ...
    xDim_across, xDim_within, min_tau, max_tau, min_eps, max_eps, ...
    min_delay, max_delay, snr, centered)
% [Ys, Xs, params] = simdata_dlag(...)
%
% Description: Generate synthetic data according to the DLAG generative 
%              model.
%
% Arguments:
%
%     N         -- int; number of sequences
%     T         -- int; number of samples per sequence
%     samplePeriod -- float; sample period (in units of time). 
%                     Assume uniform sampling.
%     yDims     -- (1 x numGroups) array; List of dimensionalities of
%                  observed data, [Y1, Y2,...]
%     xDim_across -- int; Number of across-group latents
%     xDim_within -- (1 x numGroups) array; Number of within-group latents
%                    for each group. 
%     min_tau  -- float; lower-bound of GP timescales
%     max_tau  -- float; upper-bound of GP timescales
%     min_eps  -- float; lower-bound of GP noise variances
%     max_eps  -- float; upper-bound of GP noise variances
%     min_delay -- int; lower-bound of GP timescales, in samples
%     max_delay -- int; upper-bound of GP timescales, in samples
%     snr      -- (1 x numGroups) array; List of signal-to-noise ratios,
%                 defined as trace(CC') / trace(R)
%     centered -- logical; True if data is to be zero-centered 
%                 (i.e., d = 0)
%
% Outputs:
%     Ys     -- (1 x numGroups) cell array; list of data matrices 
%               {(y1Dim x T x N), (y2Dim x T x N), ...}
%     Xs     -- (1 x numGroups) cell array; list of latent data matrices
%               {(x1Dim x T x N), (x2Dim x T x N), ...}
%     params -- structure with randomly generated observation model 
%               parameters:
%               Cs -- (1 x numGroups) cell array; List of factor loadings 
%                     {(y1Dim x x1Dim), (y2Dim x x2Dim), ...}
%               Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%                     matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%               ds -- (1 x numGroups) cell array; List of data means 
%                     {(y1Dim x 1), (y2Dim x 1), ...}
%               taus_across  -- (1 x xDim_across) array; across-group GP 
%                               timescales, in units of time 
%                               (same units as time_step)
%               taus_within -- (1 x numGroups) cell array; taus_within{i} 
%                              contains within-group GP timescales for 
%                              group i. Entries are empty for groups with 
%                              0 within-group latents.
%               eps_across -- (1 x xDim_across) array; across-group GP 
%                             noise variances
%               eps_within -- (1 x numGroups) cell array; eps_within{i} 
%                             contains within-group GP noise variances for 
%                             group i. Entries are empty for groups with 
%                             0 within-group latents.
%               delays -- (1 x numGroups) cell array; 
%                         delays{i} -- (1 x xDim_across) array; Delays 
%                         between across-group latents and group{i}.
% 
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality

    numGroups = length(yDims);
    TIME_STEP = 1; % Full temporal resolution when generating GPs
    
    % Generate GP parameters
    
    % Across-group timescales and noise variances
    taus_across = min_tau + (max_tau-min_tau).*rand(1,xDim_across); 
    % Deal with noise variances on a log scale
    eps_across = exp(log(min_eps) + (log(max_eps)- log(min_eps)).*(rand(1,xDim_across)));
    
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
    
    % Delays
    delays = cell(1,numGroups);
    for groupIdx = 1:numGroups
        if groupIdx <= 1
            % Set delays to first group to 0 time steps
            delays{groupIdx} = zeros(1,xDim_across); 
        else
            % All other groups have non-zero delays
            delays{groupIdx} = randi([min_delay, max_delay], 1, xDim_across);
        end
    end
    all_delays = cat(2,delays{:});
    largest_delay = max(abs(all_delays)); % Take the largest empirical delay
    
    % Define the sequence length of generated data *before* downsampling
    T_fullres = T * samplePeriod + 2*largest_delay;
    
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
    
    % Generate N sequences of length T according to the DLAG model
    % parameters generated above.
    Ys = cell(1, numGroups);
    Xs = cell(1, numGroups);
    for n = 1:N
        % Generate a set of latent sequences for this trial
        % Across
        Xn_across = sim_gp_latents(xDim_across, T_fullres, TIME_STEP, taus_across, eps_across);
        for groupIdx = 1:numGroups
            % If xDim_within(groupIdx) is 0, don't create within-group
            % latents
            if xDim_within(groupIdx) <= 0
                Xn_within = [];
                delay_within = [];
            else  % Otherwise, create within-group latents
                Xn_within = sim_gp_latents(xDim_within(groupIdx), T_fullres, TIME_STEP, taus_within{groupIdx}, eps_within{groupIdx});
                % We can use one function, 'generate_obs_delayed' for within-
                % and across-group latents if we just set delays for
                % within-group latents to zero.
                delay_within = zeros(1,xDim_within(groupIdx));
            end
            % Combine across- and within-group latents
            Xn = [Xn_across; Xn_within];
            % Pass these latents through pCCA model parameters to generate observed data
            [Yn, Xn_clipped] = generate_obs_delayed(Xn, Cs(groupIdx), Rs(groupIdx), ds(groupIdx), largest_delay, {[delays{groupIdx} delay_within]});
            Ys{groupIdx}(:,:,n) = Yn{1};         % Collect output
            Xs{groupIdx}(:,:,n) = Xn_clipped{1}; 
        end
    end
    
    % Collect output parameters
    params.Cs = Cs;
    params.Rs = Rs;
    params.ds = ds;
    params.taus_across = taus_across;
    params.taus_within = taus_within;
    params.eps_across = eps_across;
    params.eps_within = eps_within;
    params.delays = delays;
    
    % Downsample so that we have fewer time points
    Ys = downsample(Ys, samplePeriod);
    Xs = downsample(Xs, samplePeriod);
    
end