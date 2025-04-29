function [sa, sw] = gpspec(params, binWidth, varargin)
%
% [sa, sw] = gpspec(params, binWidth, ...)
%
% Description: Compute and (optionally) plot GP spectra over a range of 
%              frequencies specified by maxfreq.
%
% Arguments:
%
%     Required:
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
%     maxfreq    -- float; maximum frequency (in 1/time, where the units
%                   match that of binWidth) to consider when computing 
%                   cross spectra. (default: 0.5/binWidth, i.e., the
%                   maximum freuency without aliasing)
%     stepres    -- float; resolution of the computed cross spectra, in Hz 
%                   (default: 0.1).
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
%     sa   -- (numGroups x numGroups) cell array; sa{i,j} 
%             is a (1 x xDim_across) cell array, and sa{i,j}{r} is the 
%             complex-valued spectral density function between 
%             group i and group j, given by latent r. (units: magnitude 
%             would be in variance per Hz). sa is an anonymous function.
%     sw   -- (1 x numGroups) cell array; sw{i} is a (1 x xDim_within{i}) 
%             cell array, and sw{i}{r} is the within-group spectral
%             density function for group i, given by latent r. 
%             (units: variance per Hz). sw is an anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Dec 2022 -- Initial full revision.
%     18 Feb 2023 -- Added spectral Gaussian compatibility.
%     20 Feb 2023 -- Overhauled with anonymous functions and ability
%                    to compute/plot many versions of the complex-valued
%                    spectra.

showPlot.amplitude = true;
showPlot.phasedelay = true;
showPlot.groupdelay = false;
showPlot.cospec = false;
showPlot.quadspec = false;
showPlot.sqcoherency = false;
showPlot.tfgain = false;
showPlot.intspec = false;
maxfreq = 0.5/binWidth;
stepres = 0.1;
pdUnits = 'deg';
gdUnits = 'ms';
binToSec = 1/1000;
normalize = true;
assignopts(who,varargin);

binToHz = 1/binToSec; % Convert frequency in binWidth units to Hz
maxfreq = maxfreq * binToHz; % Convert maxfreq to Hz

numGroups = length(params.yDims);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;
freqSteps = (-maxfreq:stepres:maxfreq);  % Hz
freqSteps_half = (0:stepres:maxfreq);    % Hz

% Convert GP params to same units as binWidth
gp_params = getGPparams_dlag(params,binWidth);
groupParams = partitionParams_dlag(params);

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

% Compute spectral densities
sa = cell(numGroups,numGroups);   % Across-area amplitude functions
sw = cell(1,numGroups);           % Within-area amplitude functions
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        sa{groupIdx1,groupIdx2} = cell(1,xDim_across);
        for xIdx = 1:xDim_across
            D = (gp_params.DelayMatrix(groupIdx2,xIdx) - gp_params.DelayMatrix(groupIdx1,xIdx)).*binToSec; % Convert to sec
            tau = gp_params.tau_across(xIdx).*binToSec; % Convert to sec
            switch params.covType
                case 'rbf'
                    sa{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f) sqrt(2*pi*tau^2) ...
                        .*exp(-0.5*tau^2.*(2*pi.*f).^2) ...
                        .*exp(1i*2*pi*D.*f);
                case 'sg'
                    nu = gp_params.nu_across(xIdx).*binToHz; % Convert to Hz
                    sa{groupIdx1,groupIdx2}{xIdx} ...
                        = @(f)sqrt(0.5*pi*tau^2) ...
                        .*( exp(-0.5*tau^2.*(2*pi.*(f - nu)).^2) ...
                          + exp(-0.5*tau^2.*(2*pi.*(f + nu)).^2)) ...
                        .*exp(1i*2*pi*D.*f);
            end
            if ~normalize
                sa{groupIdx1,groupIdx2}{xIdx} = @(f) norm(groupParams{groupIdx1}.C(:,xIdx)) ...
                    .* sa{groupIdx1,groupIdx2}{xIdx}(f) .* norm(groupParams{groupIdx2}.C(:,xIdx));
            end
        end
        if groupIdx1 == groupIdx2
            sw{groupIdx1} = cell(1,xDim_within(groupIdx1));
            for xIdx = 1:xDim_within(groupIdx1)
                tau = gp_params.tau_within{groupIdx1}(xIdx).*binToSec; % Convert to sec
                switch params.covType
                    case 'rbf'
                        sw{groupIdx1}{xIdx} ...
                            = @(f) sqrt(2*pi*tau^2) ...
                            .*exp(-0.5*tau^2.*(2*pi.*f).^2);
                    case 'sg'
                        nu = gp_params.nu_within{groupIdx1}(xIdx).*binToHz; % Convert to Hz
                        sw{groupIdx1}{xIdx} ...
                            = @(f) sqrt(0.5*pi*tau^2) ...
                            .*( exp(-0.5*tau^2.*(2*pi.*(f - nu)).^2) ...
                                + exp(-0.5*tau^2.*(2*pi.*(f + nu)).^2));
                end
                if ~normalize
                    sw{groupIdx1}{xIdx} = @(f) norm(groupParams{groupIdx1}.C(:,xDim_across+xIdx)) ...
                        .* sw{groupIdx1}{xIdx}(f) .* norm(groupParams{groupIdx1}.C(:,xDim_across+xIdx));
                end
            end
        end
    end
end

% Plotting
if ~isempty(showPlot)
    % Across-group  
    if xDim_across > 0
        if showPlot.amplitude
            % Across-group amplitude functions
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        Amp_curr = abs(sa{groupIdx1,groupIdx2}{xIdx}(freqSteps));
                        subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                        hold on;
                        plot(freqSteps, Amp_curr, 'k-', 'linewidth', 1.5);
                        maxval = max(Amp_curr);
                        axis square;
                        xlabel('Frequency (Hz)');
                        ylabel('Amplitude (var/Hz)');
                        axis([0 maxfreq 0 1.05*maxval]);
                    end
                end
            end
        end

        if showPlot.phasedelay
            % Across-group phase functions
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        Phase_curr = angle(sa{groupIdx1,groupIdx2}{xIdx}(freqSteps)).*pdfactor;
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
            % Across-group group delay functions
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        groupdelay = @(f) gradient(angle(sa{groupIdx1,groupIdx2}{xIdx}(f)),stepres).*gdfactor; % Convert to same units as binWidth
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
            % Across-group co-spectra
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        cospec_curr = real(sa{groupIdx1,groupIdx2}{xIdx}(freqSteps));
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
            % Across-group quadrature spectra
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        quadspec_curr = imag(sa{groupIdx1,groupIdx2}{xIdx}(freqSteps));
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
            % Across-group coherency spectra
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        sqcoher = @(f) abs(sa{groupIdx1,groupIdx2}{xIdx}(f)).^2 ...
                            ./ (sa{groupIdx1,groupIdx1}{xIdx}(f) ...
                                .* sa{groupIdx2,groupIdx2}{xIdx}(f));
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
            % Across-group gains (amplitude of transfer functions)
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        sqcoher = @(f) abs(sa{groupIdx1,groupIdx2}{xIdx}(f)).^2 ...
                            ./ (sa{groupIdx1,groupIdx1}{xIdx}(f) ...
                                .* sa{groupIdx2,groupIdx2}{xIdx}(f));
                        tfgain = @(f) sqrt((sqcoher(f) ...
                            .* sa{groupIdx2,groupIdx2}{xIdx}(f)) ...
                            ./ sa{groupIdx1,groupIdx1}{xIdx}(f));
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
            % Across-group integrated spectra
            for xIdx = 1:xDim_across
                figure;
                for groupIdx1 = 1:numGroups
                    for groupIdx2 = 1:numGroups
                        cospec = @(f) real(sa{groupIdx1,groupIdx2}{xIdx}(f));
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

    end

    % Within-group
    for groupIdx = 1:numGroups
        if xDim_within(groupIdx) > 0
            if showPlot.amplitude
                % Within-group amplitude functions
                figure;
                for xIdx = 1:xDim_within(groupIdx)
                    Amp_curr = sw{groupIdx}{xIdx}(freqSteps);
                    subplot(1,xDim_within(groupIdx),xIdx);
                    hold on;
                    plot(freqSteps, Amp_curr, 'k-', 'linewidth', 1.5);
                    maxval = max(Amp_curr);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude (var/Hz)');
                    axis([0 maxfreq 0 1.05*maxval]);
                end
            end

            if showPlot.intspec
                figure;
                % Within-group integrated spectra
                for xIdx = 1:xDim_within(groupIdx)
                    subplot(1,xDim_within(groupIdx),xIdx);
                    hold on;
                    intspec = @(f) cumtrapz(f,2*sw{groupIdx}{xIdx}(f));
                    intspec_curr = intspec(freqSteps_half);
                    plot(freqSteps_half, intspec_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Frequency (Hz)');
                    ylabel('Integrated spectrum (var)');
                    axis([0 maxfreq min(intspec_curr) max([max(intspec_curr) 1.05])]);
                end
            end

        end
    end
end