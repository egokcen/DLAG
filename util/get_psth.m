function psth = get_psth(seq, varargin)
%
% psth = get_psth(seq,...)
%
% Description: Trial-average binned observations and return the 
%              peristimulus time histogram (PSTH). 
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- neural data
%
%     Optional:
%
%     spec   -- string; field name in seq to average over, if different
%               from 'y' (default: 'y')
%
% Outputs:
% 
%     psth    -- (yDim x max(T)) trial averaged activity
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision. 

spec = 'y';
assignopts(who,varargin);

Tall = [seq.T];
Tmax = max(Tall);
psth = zeros(size(seq(1).(spec),1),Tmax);
for n = 1:length(seq)    
    psth(:,1:Tall(n)) = psth(:,1:Tall(n)) + seq(n).(spec);
end

psth = psth/length(seq);

