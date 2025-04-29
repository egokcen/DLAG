function S = make_S_dlag(params, f)
%
% S = make_S_dlag(params, f)
%
% Constructs GP spectral density matrix for a given frequency.
%
% Arguments:
%
%     params -- DLAG model parameters
%     f      -- float; frequency, in units of (1/timeSteps)
%
% Outputs:
%
%     S      -- GP spectral density matrix with dimensions (xDim x xDim).
%               The matrix is diagonal.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     06 Jun 2023 -- Initial full revision.
%     13 Aug 2024 -- Only return a vector of diagonal elements, not
%                    a full diagonal matrix, plus some runtime
%                    optimizations.

numGroups = length(params.xDim_within);
xDim = params.xDim_across + sum(params.xDim_within);
% We'll use this object to track the next starting point to fill in S
xDim_totals = cumsum([params.xDim_across params.xDim_within(1:end-1)]);
% The diagonal elements of S
S = nan(xDim, 1);

% Fill in S_across
for j = 1:params.xDim_across
    switch params.covType
        case 'rbf'
            S(j) = (1 - params.eps_across(j)) .* sqrt(2*pi/params.gamma_across(j)) ...
                   .*exp(-0.5*((2*pi.*f).^2)./params.gamma_across(j)) ...
                   + params.eps_across(j);

        case 'sg'
            S(j) = (1 - params.eps_across(j)) .* sqrt(0.5*pi/params.gamma_across(j)) ...
                   .*( exp(-0.5*((2*pi.*(f - params.nu_across(j))).^2)./params.gamma_across(j)) ...
                   + exp(-0.5*((2*pi.*(f + params.nu_across(j))).^2)./params.gamma_across(j)) ) ...
                   + params.eps_across(j);
    end     
end

% Fill in S_within
for groupIdx = 1:numGroups
    for j = 1:params.xDim_within(groupIdx)
        switch params.covType
            case 'rbf'
                S(xDim_totals(groupIdx)+j,1) = (1 - params.eps_within{groupIdx}(j)) ...
                    .* sqrt(2*pi/params.gamma_within{groupIdx}(j)) ...
                    .*exp(-0.5*((2*pi.*f).^2)./params.gamma_within{groupIdx}(j)) ...
                    + params.eps_within{groupIdx}(j);

            case 'sg'
                S(xDim_totals(groupIdx)+j,1) = (1 - params.eps_within{groupIdx}(j)) ...
                    .* sqrt(0.5*pi/params.gamma_within{groupIdx}(j)) ...
                    .*( exp(-0.5*((2*pi.*(f - params.nu_within{groupIdx}(j))).^2)./params.gamma_within{groupIdx}(j)) ...
                    + exp(-0.5*((2*pi.*(f + params.nu_within{groupIdx}(j))).^2)./params.gamma_within{groupIdx}(j)) ) ...
                    + params.eps_within{groupIdx}(j);
        end   
    end
end
