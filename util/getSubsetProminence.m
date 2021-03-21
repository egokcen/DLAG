function outprom = getSubsetProminence(inprom, dimGroups_across, dimGroups_within)
%
% outprom = getSubsetProminence(inprom, dimGroups_across, dimGroups_within)
%
% Description: Return a prominence structure with only entries
%              corresponding to the specified latent groups.
%
% Arguments:
%
%     Required:
%
%     inprom -- structure containing the following fields:
%               across -- structure with across-group prominence. Has fields
%                 LL.raw    -- (1 x numDimGroups) array; prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numDimGroups) array;prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by normalized decrease in variance 
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1 x numDimGroups) cell array; each element 
%                              contains an array of across-group latent 
%                              dimensions whose joint prominence was 
%                              evaluated.
%             within -- structure with within-group prominence. Has fields
%                 LL.raw    -- (1 x numGroups) cell array; prominence of 
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numGroups) cell array; prominence of
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by normalized decrease in variance
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1 x numGroups) cell array; each element 
%                              contains a cell array of the same format as
%                              across.dimGroups, corresponding to groups of 
%                              within-group dimensions whose joint 
%                              prominence was evaluated.
%     dimGroups_across -- (1 x numKeptGroups) array; index which
%                         across-group entries in inprom to keep.
%     dimGroups_within -- (1 x numGroups) cell array; dimGroups_within{i}
%                         is a (1 x numKeptGroups) array indexing which
%                         within-group entries in inprom to keep.
%
% Outputs:
%
%     outprom -- Structure containing subset of inprom entries.
%                Same format as inprom, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 Sep 2020 -- Initial full revision.

numGroups = length(inprom.within.LL.raw); % Number of observation groups

outprom.across.LL.raw = inprom.across.LL.raw(dimGroups_across);
outprom.across.LL.upper = inprom.across.LL.upper(dimGroups_across);
outprom.across.LL.lower = inprom.across.LL.lower(dimGroups_across);
outprom.across.VE.raw = inprom.across.VE.raw(dimGroups_across);
outprom.across.VE.upper = inprom.across.VE.upper(dimGroups_across);
outprom.across.VE.lower = inprom.across.VE.lower(dimGroups_across);
outprom.across.dimGroups = inprom.across.dimGroups(dimGroups_across);

for groupIdx = 1:numGroups
    outprom.within.LL.raw{groupIdx} = inprom.within.LL.raw{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.LL.upper{groupIdx} = inprom.within.LL.upper{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.LL.lower{groupIdx} = inprom.within.LL.lower{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.VE.raw{groupIdx} = inprom.within.VE.raw{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.VE.upper{groupIdx} = inprom.within.VE.upper{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.VE.lower{groupIdx} = inprom.within.VE.lower{groupIdx}(dimGroups_within{groupIdx});
    outprom.within.dimGroups{groupIdx} = inprom.within.dimGroups{groupIdx}(dimGroups_within{groupIdx});
end
