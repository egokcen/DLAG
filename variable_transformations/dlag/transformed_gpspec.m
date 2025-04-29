function s = transformed_gpspec(V, params, binWidth, varargin)
%
% s = transformed_gpspec(V, params, binWidth, ...)
%
% Description: Given a coupled set of basis vectors, V, compute GP
%              cross spectral densities for projections that lie along 
%              those directions.
%
% Arguments:
%
%     Required:
%
%     V       -- (1 x numGroups) cell array; V{i} is a (yDims(i) x r)
%                array, where r is the number of basis vectors.
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
%                    if covType == 'sg'
%                        nu_across -- (1 x xDim_across) array; center
%                                     frequencies for spectral Gaussians;
%                                     convert to 1/time via 
%                                     nu_across./binWidth 
%                        nu_within -- (1 x numGroups) cell array; 
%                                     center frequencies for each group
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
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
%     Optional:
%
%     showPlot   -- structure containing a variety of plot options. Specify
%                   as empty ([]) to avoid plotting anything.
%                   amplitude   -- logical; plot amplitude functions
%                                  (default: true)
%                   phasedelay  -- logical; plot phase delay functions
%                                  (default: true)
%                   groupdelay  -- logical; plot group delay functions
%                                  (default: false)
%                   cospec      -- logical; plot co-spectra 
%                                  (default: false)
%                   quadspec    -- logical; plot quadrature spectra
%                                  (default: false)
%                   sqcoherency -- logical; plot squared coherency spectra
%                                  (default: false)
%                   tfgain      -- logical; plot transfer function gains
%                                  (default: false)
%                   intspec     -- logical; plot integrated (cross)-spectra
%                                  (default: false)
%                   condspec    -- logical; plot conditional spectra
%                                  (default: false)
%     maxfreq    -- float; maximum frequency (in 1/time, where the units
%                   match that of binWidth) to consider when computing 
%                   cross spectra. (default: 0.5/binWidth, i.e., the
%                   maximum freuency without aliasing)
%     stepres    -- float; resolution of the computed cross spectra, in Hz 
%                   (default: 0.1).
%     normalize  -- logical; true to normalize cross spectra; false 
%                   otherwise (default: true)
%     pdUnits    -- string; units with which to display phase functions.
%                   'deg' for degrees, 's' for time delay in seconds,
%                   'ms' for time delay in milliseconds (default: 'deg').
%     gdUnits    -- string; units with which to display group delay.
%                   's' for time delay in seconds, 'ms' for time delay in 
%                   milliseconds (default: 'ms').
%     binToSec   -- float; conversion factor from the units of binWidth
%                   to seconds (default: 1/1000, i.e., binWidth is assumed
%                   to be in ms by default)
%     normalize -- logical; set true to compute normalized spectral
%                  density functions. (default: true)
%
% Outputs:
%
%     Amp   -- (1 x r) cell array; Amp{r} is a (numGroups x numGroups) cell
%              array, and Amp{r}{i,j} is the amplitude function between
%              group i and group j. (units: variance per Hz)
%     Phase -- (1 x r) cell array; Phase{r} is a (numGroups x numGroups) 
%              cell array, and Phase{r}{i,j} is the phase function between
%              group i and group j (units in radians). (units: given by 
%              phaseUnits)
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     29 Dec 2022 -- Initial full revision.

showPlot.amplitude = true;
showPlot.phasedelay = true;
showPlot.groupdelay = false;
showPlot.cospec = false;
showPlot.quadspec = false;
showPlot.sqcoherency = false;
showPlot.tfgain = false;
showPlot.intspec = false;
showPlot.condspec = false;
maxfreq = 0.5/binWidth;
stepres = 0.1;
pdUnits = 'deg';
gdUnits = 'ms';
binToSec = 1/1000;
normalize = true;
assignopts(who,varargin);

binToHz = 1/binToSec; % Convert frequency in binWidth units to Hz
maxfreq = maxfreq * binToHz; % Convert maxfreq to Hz

numGroups = length(V);
r = size(V{1},2);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;
freqSteps = (-maxfreq:stepres:maxfreq);  % Hz
freqSteps_half = (0:stepres:maxfreq);    % Hz

% Initialize output structure
s = cell(1,r);
for j = 1:r
    s{j} = cell(numGroups); 
end

% Unit conversions for phase delay
switch pdUnits
    case 'rad'
        pdfactor = 1;
        pdylbl = 'Phase (rad)';
    case 'deg'
        pdfactor = 180 / pi;
        pdylbl = 'Phase (deg)';
    case 's'
        pdfactor = (2*pi*freqSteps).^(-1);
        pdylbl = 'Phase delay (s)';
    case 'ms'
        pdfactor = (2*pi*freqSteps./binToHz).^(-1);
        pdylbl = 'Phase delay (ms)';
end

% Unit conversions for group delay
switch gdUnits
    case 's'
        gdfactor = (2*pi).^(-1);
        gdylbl = 'Group delay (s)';
    case 'ms'
        gdfactor = (2*pi*binToSec).^(-1);
        gdylbl = 'Group delay (ms)';
end

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);
Ca = cell(1,numGroups);
Cw = cell(1,numGroups);
R = cell(1,numGroups);
for groupIdx = 1:numGroups
    R{groupIdx} = groupParams{groupIdx}.R;
    acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, cell(1,numGroups));
    Ca{groupIdx} = acrossParams.C; 
    withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:xDim_within(groupIdx)});
    Cw{groupIdx} = withinParams.C;
end

% Compute prior spectral densities
[sa_prior, sw_prior] = gpspec(params, binWidth, 'showPlot', []);

% Compute transformed cross spectra
for j = 1:r
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            s{j}{groupIdx1,groupIdx2} = @(f) 0;
            for xIdx = 1:xDim_across
                s{j}{groupIdx1,groupIdx2} = @(f) s{j}{groupIdx1,groupIdx2}(f) ...
                    + sa_prior{groupIdx1,groupIdx2}{xIdx}(f) ...
                    .*(V{groupIdx1}(:,j)'*Ca{groupIdx1}(:,xIdx)*Ca{groupIdx2}(:,xIdx)'*V{groupIdx2}(:,j));
            end
            if groupIdx1 == groupIdx2
                for xIdx = 1:xDim_within(groupIdx1)
                    s{j}{groupIdx1,groupIdx2} = @(f) s{j}{groupIdx1,groupIdx2}(f) ...
                    + sw_prior{groupIdx1}{xIdx}(f) ...
                    .*(V{groupIdx1}(:,j)'*Cw{groupIdx1}(:,xIdx)*Cw{groupIdx1}(:,xIdx)'*V{groupIdx1}(:,j)); 
                end
            end
        end
    end
end

if normalize
    
    % Compute GP covariance functions at zero lag, for normalization
    k = transformed_gpcov(V, params, binWidth, ...
                          'showPlot', false, ...
                          'normalize', false);
    % Normalize spectral densities using the GP covariances at zero lag
    for j = 1:r
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                s{j}{groupIdx1,groupIdx2} ...
                    = @(f) k{j}{groupIdx1,groupIdx1}(0).^(-0.5) ...
                      .* s{j}{groupIdx1,groupIdx2}(f) ...
                      .* k{j}{groupIdx2,groupIdx2}(0).^(-0.5);
            end
        end
    end
    
end

if ~isempty(showPlot)
    
    if showPlot.amplitude
        % Amplitude functions
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    Amp_curr = abs(s{j}{groupIdx1,groupIdx2}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, Amp_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(Amp_curr);
                    if maxval == 0
                       maxval = 1; 
                    end
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude (var/Hz)');
                    axis([0 maxfreq 0 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.phasedelay
        % Phase functions
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    Phase_curr = angle(s{j}{groupIdx1,groupIdx2}(freqSteps)).*pdfactor;
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, Phase_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(Phase_curr));
                    if maxval == 0
                       maxval = 1; 
                    end
                    line([0 maxfreq], [0 0], 'color', 'k', 'linestyle', '--');
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(pdylbl);
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.groupdelay
        % Group delay functions
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    groupdelay = @(f) gradient(angle(s{j}{groupIdx1,groupIdx2}(f)),stepres).*gdfactor; % Convert to same units as binWidth
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, groupdelay(freqSteps), 'k-', 'linewidth', 1.5);
                    maxval = max(abs(groupdelay(freqSteps)));
                    if maxval == 0
                       maxval = 1; 
                    end
                    line([0 maxfreq], [0 0], 'color', 'k', 'linestyle', '--');
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(gdylbl);
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.cospec
        % Co-spectra
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    cospec_curr = real(s{j}{groupIdx1,groupIdx2}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, cospec_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(cospec_curr));
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Co-spectrum (var/Hz)');
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.quadspec
        % Quadrature spectra
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    quadspec_curr = imag(s{j}{groupIdx1,groupIdx2}(freqSteps));
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, quadspec_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(abs(quadspec_curr));
                    if maxval == 0
                       maxval = 1; 
                    end
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Quadrature spectrum (var/Hz)');
                    axis([0 maxfreq -1.05*maxval 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.sqcoherency
        % Squared coherency spectra
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    sqcoher = @(f) abs(s{j}{groupIdx1,groupIdx2}(f)).^2 ...
                        ./ (s{j}{groupIdx1,groupIdx1}(f) ...
                            .* s{j}{groupIdx2,groupIdx2}(f));
                    sqcoher_curr = sqcoher(freqSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, sqcoher_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Sq. coherency');
                    axis([0 maxfreq 0 1.05]);
                end
            end
        end
    end
    
    if showPlot.tfgain
        % Gains (amplitude of transfer functions)
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    sqcoher = @(f) abs(s{j}{groupIdx1,groupIdx2}(f)).^2 ...
                        ./ (s{j}{groupIdx1,groupIdx1}(f) ...
                            .* s{j}{groupIdx2,groupIdx2}(f));
                    tfgain = @(f) sqrt((sqcoher(f) ...
                        .* s{j}{groupIdx2,groupIdx2}(f)) ...
                        ./ s{j}{groupIdx1,groupIdx1}(f));
                    tfgain_curr = tfgain(freqSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, tfgain_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(tfgain_curr);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(sprintf('Gain, group %d to %d', groupIdx1, groupIdx2));
                    axis([0 maxfreq 0 1.05*maxval]);
                end
            end
        end
    end
    
    if showPlot.intspec
        % Integrated spectra
        for j = 1:r
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    cospec = @(f) real(s{j}{groupIdx1,groupIdx2}(f));
                    intspec = @(f) cumtrapz(f,2*cospec(f));
                    intspec_curr = intspec(freqSteps_half);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps_half, intspec_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Integrated spectrum (var)');
                    axis([0 maxfreq min(intspec_curr) max([max(intspec_curr) 1.05])]);
                end
            end
        end
    end
    
    if showPlot.condspec
        % Conditional spectra
        for j = 1:r
            figure;
            maxvals = nan(numGroups);
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    sqcoher = @(f) abs(s{j}{groupIdx1,groupIdx2}(f)).^2 ...
                        ./ (s{j}{groupIdx1,groupIdx1}(f) ...
                            .* s{j}{groupIdx2,groupIdx2}(f));
                    condspec = @(f) sqcoher(f) .* s{j}{groupIdx2,groupIdx2}(f);
                    condspec_curr = condspec(freqSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    plot(freqSteps, condspec_curr, 'k-', 'linewidth', 1.5);
                    maxvals(groupIdx1,groupIdx2) = max(condspec_curr);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel(sprintf('Amplitude, group %d due to %d', groupIdx2, groupIdx1));
                end
            end
            % Set the vertical axis scales equal across rows
            maxval = max(maxvals,[],2);
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    axis([0 maxfreq 0 1.05*maxval(groupIdx1)]);
                end
            end
        end
    end
    
end