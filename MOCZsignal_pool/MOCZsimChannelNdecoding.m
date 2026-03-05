function [P_rec] = MOCZsimChannelNdecoding(simParams, SNR, x_pb, P)
%MOCZSIMCHANNELNDECODING Summary of this function goes here
%   Detailed explanation goes here
% assuming x_pb is a vector that is on pass band - contain zadoffchu and
% dataload
    K = simParams.K;
    B = simParams.B;
    P_rec = zeros(K, B);

    %% 1. Adding Channel
    Tsym = simParams.Tsym;
    Fs = simParams.Fs;
    L = simParams.L;
    
    %option 1 - Rayleigh Channel
%     pathDelays = (0:L-1) * Tsym; % L taps spaced by 1 symbol each
%     pathGains = [0 -15 -30 -45 -60];        % Average path gains in dB
%     
%     fadingChannel = comm.RayleighChannel(...
%         'SampleRate', Fs, ...
%         'PathDelays', pathDelays, ...
%         'AveragePathGains', pathGains, ...
%         'MaximumDopplerShift', 0); % Adjust Doppler for motion/fading speed
% %     Low Doppler (e.g., 0 to 2 Hz): Represents walking speed or a mostly stationary environment.
% %     The fading channel stays relatively stable from symbol to symbol.
% 
%     % The output will be longer than the input due to the delay spread -
%     % but for some reason, its not..
%     chanOutput = fadingChannel(x_pb);
%     
%     % Note: To see the "extra" symbols, you may need to flush the channel 
%     % or append zeros to the input to allow the multi-path tails to exit.
%     release(fadingChannel);

    %option 2 - static fir
    h_channel = [0.75, -0.35, 0.15, -0.07, 0.003];
    sps = Tsym * Fs; %samples per symbol
    h_channel_upsmp = upsample(h_channel, sps);
    h_channel_upsmp = h_channel_upsmp / norm(h_channel_upsmp); %normalize
    chanOutput = filter(h_channel_upsmp, 1, x_pb);
%     chanOutput = x_pb;
    
    x_pb_noisy = awgn(chanOutput, SNR ,'measured');
    
    %% 2. Decoding
    N = simParams.zadoffChuPair(2);
    R = simParams.R;
    theta_c = simParams.theta_c;

    [~, x_decoded] = pulseSHP(x_pb_noisy, simParams, 'demodulate');
    % x_decoded is a matrix with B colls, and N+((K+1)+L)*B rows 

    % after real measurments, here there is a cross-correlation with
    % the zadoffchu sequence in order to sync. can be some func:
    % idx = findZC(x_decoded, simParams.zadoffChuPair) - idx will be the
    % idx of 1st element after the zadOffChu (1st bit of guard)
    % x_decoded = x_decoded(idx + L :end) (L is guard)
    % In simulation though we know the 1st simParams.N bits are the zadoffchu
    % sequence:
    x_decoded = x_decoded(N + L +1 : end);
    % now, x_decoded is a corrupted version of [B sequences of
    % (K+1 data bits, L guard bits)]. must check that its truly has a total
    % of B*(K+1+L) elements!
    x_decoded_matrix = reshape(x_decoded, K + 1 + L, B);
    % ommiting the guard bits:
    x_decoded_matrix = x_decoded_matrix(1 : K + 1, :);

    for b = 1 : B
        x_vec_decoded = x_decoded_matrix(:, b);
        message_rec = DiZeT(R, theta_c, x_vec_decoded, K);
        P_rec(:, b) = message_rec;
        %if want visualization:
        if simParams.figFlag
            decoded_alphas = roots(flip(x_vec_decoded'));
            % Get current message bits and extract oridginal alphas:
            M = P(:, b);
            alphas = ((1-M)*(R^(-1)) + M * R) .* exp(1i*theta_c);

            figure;
            % Plot decoded zeros - save handle 'h1'
             h1 = scatter(real(decoded_alphas), imag(decoded_alphas), 50, 'filled', 'b');
            grid on; axis equal; hold on;
            % Plot encoded zeros - save handle 'h2'
            h2 = scatter(real(alphas), imag(alphas), 50, 'filled', 'g');
            % Draw circles (we don't save handles because we don't want them in the legend)
            thetas_for_plot = linspace(0, 2*pi, 300);
            plot(R * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
            plot((1/R) * exp(1i*thetas_for_plot), 'r', 'LineWidth', 1);
            plot(exp(1i*thetas_for_plot), 'k--', 'LineWidth', 0.5);
            % Explicitly call legend only for h1 and h2
            legend([h1, h2], {'Decoded Zeros', 'Original Zeros'});
            title('Huffman BMOCZ encoded and decoded zeros');
            xlabel('Re(z)');
            ylabel('Im(z)');
            hold off;
        end
    end
end

