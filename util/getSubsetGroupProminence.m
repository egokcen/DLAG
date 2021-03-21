function outprom = getSubsetGroupProminence(inprom, dimGroups_across, dimGroups_within)
%
% outprom = getSubsetGroupProminence(inprom, dimGroups_across, dimGroups_within)
%
% Description: Return a group prominence structure with only entries
%              corresponding to the specified latent groups.
%
% Arguments:
%
%     Required:
%
%     inprom -- (1 x numGroups) cell array; inprom{i} contains a structure
%               with the bootstrapped prominence of each group of latent 
%               variables for observation group i. See bootstrapProminence 
%               for details on the format of these structures.
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

numGroups = length(inprom); % Number of observation groups

% Initialize output structure
outprom = inprom;

for groupIdx = 1:numGroups
    outprom{groupIdx}.across.LL.raw = inprom{groupIdx}.across.LL.raw(dimGroups_across);
    outprom{groupIdx}.across.LL.upper = inprom{groupIdx}.across.LL.upper(dimGroups_across);
    outprom{groupIdx}.across.LL.lower = inprom{groupIdx}.across.LL.lower(dimGroups_across);
    outprom{groupIdx}.across.VE.raw = inprom{groupIdx}.across.VE.raw(dimGroups_across);
    outprom{groupIdx}.across.VE.upper = inprom{groupIdx}.across.VE.upper(dimGroups_across);
    outprom{groupIdx}.across.VE.lower = inprom{groupIdx}.across.VE.lower(dimGroups_across);
    outprom{groupIdx}.across.dimGroups = inprom{groupIdx}.across.dimGroups(dimGroups_across);

    outprom{groupIdx}.within.LL.raw{1} = inprom{groupIdx}.within.LL.raw{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.LL.upper{1} = inprom{groupIdx}.within.LL.upper{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.LL.lower{1} = inprom{groupIdx}.within.LL.lower{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.VE.raw{1} = inprom{groupIdx}.within.VE.raw{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.VE.upper{1} = inprom{groupIdx}.within.VE.upper{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.VE.lower{1} = inprom{groupIdx}.within.VE.lower{1}(dimGroups_within{groupIdx});
    outprom{groupIdx}.within.dimGroups{1} = inprom{groupIdx}.within.dimGroups{1}(dimGroups_within{groupIdx});
end
