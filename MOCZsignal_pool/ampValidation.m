% validation script - CCPF of PAPR graph
clear;
%% Part A - obtaining values of PAPR of all possible permutations

% instead of this, can load simParams:
%%%%%
K = 12;
lambda = 1;
Kidxs = 0:K-1;
R = sqrt(1+2*lambda*sin(pi/K));
theta_c = ((2*pi) * (Kidxs / K))';
%%%%%

bits = ff2n(K);
PAPR = zeros(2^K,1);


for ii = 1 : 2^K
    % Calculate the PAPR for the current bit sequence:
    currentBits = bits(ii, :);
    M = currentBits';
    % 1. Generate Alphas (Roots)
    alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);
    % 2. Generate Coefficients
    x_un = flip(poly(alphas))'; %finding K+1 coefficients and flip - now the free coefficent is first
    x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1). now total energy of signal will be K+1

    PAPR(ii) = 10 * log10(max(abs(x).^2) / mean(abs(x).^2));
end

%% Part B - displaying CCPF of PAPR

% CCDF = 1 - histogram(PAPR, 'Normalization', 'cdf');

papr_sorted = sort(PAPR); % Sort values in ascending order
N = length(PAPR);

% Calculate the probability P(PAPR > threshold)
% For the i-th sorted value, there are (N - i) values greater than it.
prob = (N:-1:1) / N; 

% 3. Plot using a log scale for the Y-axis
figure;
semilogy(papr_sorted, prob, 'LineWidth', 2);
grid on;
xlabel('PAPR_0 [dB]');
ylabel('Prob(PAPR > PAPR_0)');
title('CCDF of PAPR');