function tau = log_tau_range(n_tau,tau_range)
tau = 1:n_tau;
tau = exp(tau/n_tau);
tau = tau - tau(1);
tau = tau / tau(end);
tau = tau * (tau_range(2) - tau_range(1)) + tau_range(1);
end

