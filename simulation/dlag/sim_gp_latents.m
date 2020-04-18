function X = sim_gp_latents(xDim, T, time_step, taus, eps)
% 
% X = sim_gp_latents(xDim, T, time_step, taus, eps)
%
% Description: Generate xDim independent sequences of length T samples,
%              according to a zero-mean Gaussian Process with squared 
%              exponential kernel.
%
% Arguments:
%     xDim      -- int; number of latents sequences to generate
%     T         -- int; number of samples
%     time_step -- float; sample period (in units of time). 
%                  Assume uniform sampling.
%     taus      -- (1 x xDim) array; timescales for each latent gaussian 
%                  process (assuming RBF kernel)
%     eps       -- (1 x xDim) array; noise variances of each gaussian
%                  process (assuming normalized RBF kernel, so that k(0)=1)
%
% Outputs:
%     X -- (xDim x T) array; simulated latent sequences
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.
    
    X = zeros(xDim, T);
    for i = 1:xDim
        Ki = rbf(T, taus(i), time_step, eps(i));
        X(i,:) = mvnrnd(zeros(1,T), Ki);
    end

end