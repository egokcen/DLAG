function K = rbf(T, tau, time_step, eps)
% Helper function to compute a radial basis function (RBF) or squared exponential kernel for
% a specified number of time points.

    gamma = (time_step / tau)^2; % The 'timescale' is now in terms of sample indices
    Ts = (1:T)';
    Tdiff = squareform(pdist(Ts, 'squaredeuclidean'));
    K = (1 - eps) .* exp(-(gamma/2) .* Tdiff) + eps .* eye(T);
end