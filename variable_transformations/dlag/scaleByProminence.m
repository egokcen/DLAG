function seqScaled = scaleByProminence(seq, prom, xDim_across, xDim_within, varargin)
%
% seqScaled = scaleByProminence(seq, prom, xDim_across, xDim_within...)
%
% Description: Rescale inferred latent trajectories according to their
%              prominence within each group. This rescaling is for
%              visualization purposes only.
%
%              NOTE: This function assumes that prominence was evaluated on
%                    each latent variable individually, not on groups of
%                    latents jointly.
%
% Arguments:
%
%     Required:
%
%     seq -- data structure whose i-th entry (corresponding to the i-th
%            trial) has the following (relevant) field:
%            xsm -- ((numGroups*xDim) x T) array; posterior mean at each 
%                   timepoint
%
%     prom -- (1 x numGroups) cell array; prom{i} contains a structure with
%             the bootstrapped prominence of each group of latent variables
%             for observation group i. See bootstrapProminence for details
%             on the format of these structures.
%
%     xDim_across -- int; Number of across-group dimensions.
%     
%     xDim_within -- (1 x numGroups) array; Number of within-group
%                    dimensions for each group.
%
%     Optional:
%
%     metric -- string; Specify which measure of prominence use: 'LL'
%               or 'VE' (default: 'VE')
%
% Outputs:
%
%     seqScaled -- rescaled version of seq with the following changed field:
%                  xsm -- ((numGroups*xDim) x T) array; rescaled posterior 
%                         mean at each timepoint
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Sep 2020 - Initial full revision.

% Optional arguments
metric = 'VE';
assignopts(who,varargin);

N = length(seq); % Number of trials
numGroups = length(xDim_within); % Number of observation groups

% Initialize output structure
seqScaled = seq;

% Partition latents according to observation group
groupSeq = partitionObs(seq, xDim_across + xDim_within, 'datafield', 'xsm');

% Partition latents into across- and within-group sets
groupSeq_across = cell(1,numGroups);
groupSeq_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    [seqAcross, seqWithin] = partitionLatents_meanOnly(groupSeq{groupIdx}, xDim_across, xDim_within(groupIdx));
    groupSeq_across{groupIdx} = seqAcross;
    groupSeq_within{groupIdx} = seqWithin{1};
end

% Turn prominence values into normalized scale factors
scaleFactors_across = cell(1,numGroups);
scaleFactors_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Combine individual-group prominences of across- and within-group latents
    prom_all = [prom{groupIdx}.across.(metric).raw prom{groupIdx}.within.(metric).raw{1}];
    % Convert across- and within-group prominences into normalized scale
    % factors
    scaleFactors_across{groupIdx} = prom{groupIdx}.across.(metric).raw ./ max(prom_all);
    scaleFactors_within{groupIdx} = prom{groupIdx}.within.(metric).raw{1} ./ max(prom_all);
end

% Scale latent trajectories by the normalized scale factors, and place them
% in the output structure.
for n = 1:N
    seqScaled(n).xsm = [];
    for groupIdx = 1:numGroups
        if xDim_across > 0
            seqScaled(n).xsm = [seqScaled(n).xsm; groupSeq_across{groupIdx}(n).xsm .* scaleFactors_across{groupIdx}'];
        end
        if xDim_within(groupIdx) > 0
            seqScaled(n).xsm = [seqScaled(n).xsm; groupSeq_within{groupIdx}(n).xsm .* scaleFactors_within{groupIdx}'];
        end
    end    
end