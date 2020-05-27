function prom = evalProminence(seq, params, dimGroups_across, dimGroups_within)
%
% prom = evalProminence(seq, params, dimGroups_across, dimGroups_within)
%
% Description: Evaluate the "prominence" of each group of latents given in 
%              dimGroups_across and dimGroups_within. "Prominence" here
%              is defined as the relative decrease in performance that
%              results from removing a group of latents from the full
%              model.
%
% Arguments:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
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
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%     dimGroups_across -- (1 x numDimGroups) cell array; each element 
%                         contains an array of across-group latent dimensions
%                         whose joint prominence will be evaluated.
%     dimGroups_within -- (1 x numGroups) cell array; each element 
%                         contains a cell array of the same format as
%                         dimGroups_across, corresponding to groups of 
%                         within-group dimensions whose joint prominence
%                         will be evaluated.
%
% Outputs:
%
%     prom -- structure containing the following fields:
%             across -- structure with across-group prominence. Has fields
%                 LL        -- (1 x numDimGroups) array; LL(i) contains
%                              prominence of dimGroups_across(i), measured
%                              by decrease in log-likelihood relative to 
%                              the full model.
%                 VE        -- (1 x numDimGroups) array; VE(i) contains
%                              prominence of dimGroups_across(i), measured
%                              by normalized decrease in variance explained
%                              relative to the full model.
%                 dimGroups -- (1x numDimGroups) cell array; Exactly
%                              dimGroups_across. Entries correspond to
%                              above fields.
%             within -- structure with within-group prominence. Has fields
%                 LL        -- (1 x numGroups) cell array; LL{i} contains
%                              prominence of dimGroups_within{i}, measured
%                              by decrease in log-likelihood relative to 
%                              the full model.
%                 VE        -- (1 x numGroups) cell array; VE{i} contains
%                              prominence of dimGroups_within{i}, measured
%                              by normalized decrease in variance explained
%                              relative to the full model.
%                 dimGroups -- (1x numGroups) cell array; Exactly
%                              dimGroups_within. Entries correspond to
%                              above fields.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 May 2020 -- Initial full revision.

numGroups = length(dimGroups_within);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
yDims = params.yDims;

% Initialize output structure
prom.across.LL = nan(1,length(dimGroups_across));
prom.across.VE = nan(1,length(dimGroups_across));
prom.across.dimGroups = dimGroups_across;

prom.within.LL = cell(1,numGroups);
prom.within.VE = cell(1,numGroups);
for groupIdx = 1:numGroups
    prom.within.LL{groupIdx} = nan(1,length(dimGroups_within{groupIdx})); 
    prom.within.VE{groupIdx} = nan(1,length(dimGroups_within{groupIdx})); 
end
prom.within.dimGroups = dimGroups_within;

% Structures for when we want to keep all dimensions, within or across
allAcross = 1:xDim_across;
allWithin = cell(1,numGroups);
for groupIdx = 1:numGroups
    allWithin{groupIdx} = 1:xDim_within(groupIdx); 
end

% Evaluate the likelihood of the full model, to provide a baseline.
[~, LL_full] = exactInferenceWithLL_dlag(seq, params, 'getLL', true);

% Evaluate the variance explained by the full model, to provide a baseline.
[~, ~, R2_full] = denoise_dlag(seq, params);

% Evaluate across-group prominence
for dimIdx = 1:length(dimGroups_across)
    dimGroup = dimGroups_across{dimIdx};
    % Take all latent dimensions but the ones we wish to evaluate
    keptAcross = setdiff(allAcross,dimGroup);
    subparams = getSubsetParams_dlag(params, keptAcross, allWithin);
    % Workaround to handle zero across- and within-group dimensions
    if subparams.xDim_across == 0 && ~any(subparams.xDim_within)
        seq_static = seq2pcca(seq, yDims, 'datafield', 'y');
        obsBlockIdxs = get_block_idxs(yDims);
        params_static.Rs = mat2blocks(subparams.R, obsBlockIdxs);
        params_static.ds = cell(1,numGroups);
        for groupIdx = 1:numGroups
            obsBlock = obsBlockIdxs{groupIdx}(1):obsBlockIdxs{groupIdx}(2);
            params_static.ds{groupIdx} = subparams.d(obsBlock);
        end
        % LL
        rGroups = [1 2];
        [LL_sub, ~, ~] = indepGroupEval(seq_static, params_static, rGroups);
        % R^2
        seqAll = cat(1, seq_static{:});
        seqPred = repmat(subparams.d, [1 size(seqAll,2)]);
        RSS = sum( sum( ( seqAll - seqPred ).^2, 1 ) );
        TSS = sum( sum( ( seqAll - repmat( mean(seqAll,2), [1 size(seqAll,2)] ) ).^2, 1 ) );
        R2_sub = 1 - RSS / TSS;
    else
        % Evaluate prominence normally
        [~, LL_sub] = exactInferenceWithLL_dlag(seq, subparams, 'getLL', true);
        [~, ~, R2_sub] = denoise_dlag(seq, subparams);
    end
    % Store results
    prom.across.LL(dimIdx) = LL_full - LL_sub;
    prom.across.VE(dimIdx) = (R2_full - R2_sub) / R2_full;
end

% Evaluate within-group prominence
for groupIdx = 1:numGroups
    keptWithin = allWithin;
    for dimIdx = 1:length(dimGroups_within{groupIdx})
        dimGroup = dimGroups_within{groupIdx}{dimIdx};
        % Take all latent dimensions but the ones we wish to evaluate
        keptWithin{groupIdx} = setdiff(allWithin{groupIdx},dimGroup);
        subparams = getSubsetParams_dlag(params, allAcross, keptWithin);
        % Evaluate prominence
        [~, LL_sub] = exactInferenceWithLL_dlag(seq, subparams, 'getLL', true);
        [~, ~, R2_sub] = denoise_dlag(seq, subparams);
        % Store results
        prom.within.LL{groupIdx}(dimIdx) = LL_full - LL_sub;
        prom.within.VE{groupIdx}(dimIdx) = (R2_full - R2_sub) / R2_full;
    end
end