function [P, signal_pb_total, x_bb, group_delay] = generateFunc(simParams)
% -- Params --
% 1. MOCZ
B = simParams.B;
L = simParams.L;
K = simParams.K;

% 2. zadOffChu (sync):
N = simParams.zadoffChuPair(2);
u = simParams.zadoffChuPair(1);


% P is a binary packet. there are B messages (columns) and K bits in each
% message
P = randi([0,1], K, B); 
signal_pb_total = [];

% -- Huffman BMOCZ --
R = simParams.R;
theta_c = simParams.theta_c;

% for sync:
zc = zadoffChuSeq(u, N);
% Initialize an empty vector for all baseband symbols
symbols_total = [];
% Add Zadoff-Chu sequence and its guard to the beginning
zc_with_guard = [zc; zeros(L, 1)];
symbols_total = [symbols_total; zc_with_guard];

% Process each message in the packet
for b = 1:B
    % Get current message bits
    M = P(:, b);
    
    % Generate Alphas (Roots)
    alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);
    
    % Generate Coefficients
    x_un = flip(poly(alphas))'; %finding K+1 coefficients and flip - now the free coefficent is first
    x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1). now total energy of signal will be K+1
    
    x_with_guard = [x; zeros(L, 1)]; %Guard of zeros (channel of L taps)
    
    symbols_total = [symbols_total; x_with_guard];

end

[signal_pb_total, ~, group_delay] = pulseSHP(symbols_total, simParams, 'modulate');
x_bb = signal_pb_total;
end