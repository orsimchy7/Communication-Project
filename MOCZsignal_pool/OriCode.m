%ORY SCHAUL & AMIT AMRAM
%This script runs a Monte Carlo simulation comparing all decoding methods evaluated 
%in this project, including ML, DiZeT, Greedy, and the modulation schemes FSK and DPSK.
clc; clear; close all;

%%  creating MOCZ parameters
K = 7;           % Message length in bits
L = 3;           % Channel length
N = K + L ;   % Convolution length
R = sqrt(1 + sin(2*pi/K));  % Radius based on Huffman zero placement
gain = sort(rand(1,L-1), 'descend'); %maybe write gain in DB
h = [1, sqrt(gain).*exp(1i*2*pi*rand(1,L-1))]; %time domain taps. h zeros are root(h)

SNR = -10:2:30;         % SNR values in dB
numIter = 10000;          % Iterations per SNR point
BER = zeros(length(SNR),1,5);
angles = 2*pi*(0:K-1)/K;

%Greedy Variables
zeros_bank = zeros(1,2*K); %2K optional zeros of polynom
zeros_bank(1,1:K) = R*exp(1i * angles(1:K)); %ONES
zeros_bank(1,(K+1):2*K) = (1/R)*exp(1i * angles(1:K)); %ZEROS
V_bank = zeros(N, 2*K); %all Vj vectors
for col = 1:2*K
    V_bank(:, col) = conj(zeros_bank(col)).^(0:N-1);
end

%ML variables
% Precompute all 2^K binary combinations and their corresponding zero vectors
combinations = dec2bin(0:2^K-1) - '0';   %matrix of all 2^k possible messages
V_all = zeros(N,K,2^K); %all possible vandermonde matrixes

alpha_candidates = zeros(2^K, K); %all possible alphas - mapping
for idx = 1:2^K
    for k = 1:K
        if combinations(idx, k) == 1
            alpha_candidates(idx, k) = R * exp(1i * angles(k));
        else
            alpha_candidates(idx, k) = (1/R) * exp(1i * angles(k));
        end
    end
    
    alpha = alpha_candidates(idx, :);
    V = zeros(N, K);
    for col = 1:K
        V(:, col) = conj(alpha(col)).^(0:N-1);
    end

    V_all(:,:,idx) = V;

    C = V' * V;
    [U, S_diag, ~] = svd(C);
    Cinv_all(:,:,idx) = U * diag(1 ./ sqrt(diag(S_diag))) * U';
    % Cinv_all(:,:,idx) = inv(sqrtm(C));
end
%% FSK and DPSK parameters
M = 2;
Fs = 1;
nsamp = 8;
freqsep = Fs/nsamp;    % q = 1

h_up = zeros(1, (numel(h)-1)*nsamp + 1);
for ell = 0:numel(h)-1
    h_up(1 + ell*nsamp) = h(1 + ell);
end

%% Simulation
for i = 1:length(SNR)
    BER_i = zeros(numIter, 1,5);
    for j = 1:numIter
        % Transmitter 
        bits = randi([0, 1], K, 1);  % Random bit message
        zeros_tx = zeros(K,1);
        for k = 1:K
            if bits(k) == 1
                zeros_tx(k) = R * exp(1i * angles(k));
            else
                zeros_tx(k) = (1/R) * exp(1i * angles(k));
            end
        end
        poly_coeffs = fliplr(poly(zeros_tx));   
        % Channel + Noise 
        channel_output = conv(poly_coeffs, h);
        noisy_signal = awgn(channel_output, SNR(i), 'measured');
        y = noisy_signal(:);  % Ensure y is a column vector of size Nx1

        %MOCZ
        % decoding with ML
        min_cost = inf;
        bits_rx_ML = zeros(K, 1);
        
        for idx = 1:2^K
            V = V_all(:,:,idx);          % Vandermonde matrix 
            Cinv = Cinv_all(:,:,idx);    % Inverse sqrt Gram matrix (K×K)

            Vy = V' * y;                  % (K×1)
            projection = Cinv * Vy;      % (K×1)
            cost = norm(projection)^2;

            if cost < min_cost
                min_cost = cost;
                bits_rx_ML = combinations(idx, :).';
            end
        end        

        % Count bit errors
        [~,BER_i(j,:,1)] = biterr(bits,bits_rx_ML);
        
        %DiZeT
        for k = 1:K %decoding
            zero_k = R*exp(2i*pi*(k-1)/K);
            conj_zero_k = (1/R)*exp(2i*pi*(k-1)/K);
            Y_zk = abs(polyval(fliplr(noisy_signal), zero_k));
            Y_zk_conj = abs(polyval(fliplr(noisy_signal), conj_zero_k));
            if Y_zk > (R^(N))*Y_zk_conj
                bits_rx_DiZeT(k) = 0;
            else
                bits_rx_DiZeT(k) = 1;
            end
        end
        [~,BER_i(j,:,2)] = biterr(bits,bits_rx_DiZeT.');


        %Greedy algorithm
        %Initialization with DiZeT
        bits_rx_greedy = zeros(K,1);
        zero_1 = R;
        conj_zero_1 = (1/R);
        Y_zk = abs(polyval(fliplr(noisy_signal), zero_1));
        Y_zk_conj = abs(polyval(fliplr(noisy_signal), conj_zero_1));
        if Y_zk > (R^(N))*Y_zk_conj
            bits_rx_greedy(1) = 0;
            V_temp = V_bank(:,K+1);
        else
            bits_rx_greedy(1) = 1;
            V_temp = V_bank(:,1);
        end    
        mask = zeros(1,2*K);
        mask(1) = 1;
        mask(K+1) = 1;
        P_i = V_temp*inv(V_temp'*V_temp)*V_temp'; %NxN
        P_i_orth = eye(N,N)-P_i;
        
        % greedy selection loop
        for t = 2:K
            min_cost = inf;
            min_l = 0;
            for l = 1:2*K
                if mask(l) == 1
                    continue
                end
                u_i = (P_i_orth*V_bank(:,l))/norm((P_i_orth*V_bank(:,l)));
                projection = (abs(u_i'*y))^2;
                if projection < min_cost
                    min_cost = projection;
                    min_l = l;
                end
            end
            V_temp = [V_temp V_bank(:,min_l)];
            if min_l > K
                k_idx = min_l - K;
                mask(min_l) = 1;
                mask(min_l - K) = 1;
                bit_val = 0;
            else
                k_idx = min_l;
                mask(min_l) = 1;
                mask(min_l+K) = 1;
                bit_val = 1;
            end
            bits_rx_greedy(k_idx) = bit_val;

            P_i = V_temp*inv((V_temp)'*V_temp)*V_temp';    
            P_i_orth = eye(N,N)-P_i;
        end
        [~,BER_i(j,:,3)] = biterr(bits,bits_rx_greedy);
        
        % %FSK
        % txFSK = fskmod(bits, M, freqsep, nsamp, Fs);   % column vector
        % yFSK    = conv(txFSK, h_up);                        % length K*nsamp + L - 1
        % rFSK_clean = yFSK(1:K*nsamp);
        % rFSK = awgn(rFSK_clean, SNR(i), 'measured');
        % dataFSK = fskdemod(rFSK, M, freqsep, nsamp, Fs);
        % [~, BER_i(j,:,4)] = biterr(bits, dataFSK);        
        % 
        % 
        % % DPSK 
        % txDPSK = dpskmod(bits, M);
        % txDPSK_z = [txDPSK;zeros(numel(h)-1,1)];
        % rxDPSK_channel = conv(txDPSK_z,h);
        % rxDPSK = rxDPSK_channel(1:numel(txDPSK));           % trim only the zero-flush tail
        % rxDPSK = awgn(rxDPSK, SNR(i), 'measured');
        % dataDPSK = dpskdemod(rxDPSK, M);
        % [~,BER_i(j,:,5)] = biterr(bits,dataDPSK);

    end
    
    BER(i,:,1) = mean(BER_i(:,1,1));
    BER(i,:,2) = mean(BER_i(:,1,2));
    BER(i,:,3) = mean(BER_i(:,1,3));
    BER(i,:,4) = mean(BER_i(:,1,4));
    BER(i,:,5) = mean(BER_i(:,1,5));

end
%%
% --- Plotting BER vs SNR ---
figure;
semilogy(SNR, BER(:,:,1), 'r-o', 'LineWidth', 1.5,'DisplayName', 'ML');
xlabel('SNR [dB]');
ylabel('Bit Error Rate (BER)');
grid on;
hold on
semilogy(SNR, BER(:,:,2), 'b-o', 'LineWidth', 1.5,'DisplayName', 'DiZeT');
semilogy(SNR, BER(:,:,3), 'g-o', 'LineWidth', 1.5,'DisplayName', 'Greedy');
semilogy(SNR, BER(:,:,4), 'y-o', 'LineWidth', 1.5,'DisplayName', 'FSK');
semilogy(SNR, BER(:,:,5), 'm-o', 'LineWidth', 1.5,'DisplayName', 'DPSK');
title('Decoding methods comparison');
legend('show', 'Location', 'northeast');


%%
sym0 = fskmod(0, M, freqsep, nsamp, Fs);
sym1 = fskmod(1, M, freqsep, nsamp, Fs);
Nfft = 1024; f = (-Nfft/2:Nfft/2-1)/Nfft * Fs;
S0 = fftshift(abs(fft(txFSK, Nfft)));
S1 = fftshift(abs(fft(sym1, Nfft)));

figure;
plot(f, S0/max(S0), 'r', 'LineWidth',1.2); hold on
%plot(f, S1/max(S1), 'b', 'LineWidth',1.2); grid on
legend('bit=0','bit=1'); xlabel('Normalized frequency'); ylabel('Mag (norm)');
title('One-symbol spectra: two tones at \pm freqsep/2');

%%
% --- FFT settings ---
Nfft = 4096;
Fs_sym = 1;     % symbol-rate "sampling freq"
Fs_samp = 1;    % your Fs for nsamp samples/symbol (consistent with fskmod)

% --- Symbol-rate response H_sym(f) ---
Hsym = fftshift(fft(h, Nfft));
f_sym = (-Nfft/2:Nfft/2-1)/Nfft * Fs_sym;

figure; 
subplot(2,1,1);
plot(f_sym, 20*log10(abs(Hsym)+1e-12), 'LineWidth', 1.2); grid on
xlabel('Normalized frequency (symbol rate)'); ylabel('|H_{sym}(f)| [dB]');
title('Symbol-rate channel response  H_{sym}(f)');

subplot(2,1,2);
plot(f_sym, angle(Hsym), 'LineWidth', 1.0); grid on
xlabel('Normalized frequency (symbol rate)'); ylabel('\angle H_{sym}(f) [rad]');