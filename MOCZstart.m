function [M, message_rec] = MOCZstart(simParams, SNR)

    % -- Params --
    B = simParams.B;
    L = simParams.L;
    K = simParams.K + L; %add L-1 only if adding a guard of L taps
    Tsym = simParams.Tsym;
    fc = simParams.fc;
    Fs = simParams.Fs;
    lambda = simParams.lambda;
    figFlag = simParams.figFlag;
    L = simParams.L;

    Kidxs = 1 : (K);
    
    M = randi([0,1],K - (L),B);
    zeroPad = zeros(L, B);
    M = [M;zeroPad];
    % M_centered = 2*M -1;
    
    % -- Huffman BMOCZ --
    R = sqrt(1+2*lambda*sin(pi/K));
    theta_c = ((2*pi) * (Kidxs / K))';
    
    alphas = ((1-M)*(R^(-1)) + M * R).*exp(1i*theta_c);
    
    x_un = flip(poly(alphas))';
    x = (x_un / norm(x_un)) * (K+1); % normalize to K+1
    
    [x_pb, ~] = pulseSHP(x, simParams, 'modulate');

%     t = (0:1/Fs:Tsym*(K+1 +group_delay/sps) - 1/Fs)';
%     x_pb = sqrt(2) * real(x_mod.*exp(1i*2*pi*fc*t));

%     % -- Visualization --
%     if figFlag
%         figure(1);
%         subplot(2,1,1);
%         stem(h_coeff, 'filled'); title('Raised Cosine Impulse Response'); grid on;
%         subplot(2,1,2);
%         t = 1:sps*(K+1); % Time axis for 9 symbols
%         sig_segment = x_mod(t);
%         % plot3(X, Y, Z) where X is time, Y is Real, Z is Imaginary
%         plot3(t, real(sig_segment), imag(sig_segment), 'b', 'LineWidth', 1);
%         hold on;
%         % Add the discrete symbols as points in 3D space
%         t_symbols = 1:sps:sps*(K+1);
%         stem3(t_symbols, real(x), imag(x), 'r', 'LineWidth', 2);
%         view(-25, 30); % Rotate for a better 3D perspective
%         xlabel('Time (Samples)'); ylabel('Real (In-phase)'); zlabel('Imaginary (Quadrature)');
%         title('3D Complex Pulse Shaped Signal');
%         grid on;
%         hold off;
%     end

    %% adding channel
    pathDelays = (0:4) * Tsym; % 5 taps spaced by 1 symbol each
    pathGains = [0 -2 -3 -6 -10];        % Average path gains in dB
    
    % Create the Rayleigh Channel Object
    fadingChannel = comm.RayleighChannel(...
        'SampleRate', Fs, ...
        'PathDelays', pathDelays, ...
        'AveragePathGains', pathGains, ...
        'MaximumDopplerShift', L); % Adjust Doppler for motion/fading speed
    % The output will be longer than the input due to the delay spread
    chanOutput = fadingChannel(x_pb);
    
    % Note: To see the "extra" symbols, you may need to flush the channel 
    % or append zeros to the input to allow the multi-path tails to exit.
    release(fadingChannel);
    x_pb_noisy = awgn(chanOutput, SNR ,'measured');
    
    %% Decoding
    % part A - band pass filtering
    % b_bp = designfilt('bandpassfir', 'FilterOrder', 100, ...
    %     'CutoffFrequency1', fc- BW, 'CutoffFrequency2', fc + BW, 'SampleRate', Fs);
    % x_pb_filtred = filtfilt(b_bp.Coefficients, 1, x_pb_noisy);
    % freqz(b_bp.Coefficients,1,1024,Fs);
    
    %part b - shifting to 0 freq
    % mixed_x = x_pb_filtred.*exp(-1i*2*pi*fc*t);
    % 
    % %part c - LP filter
    % f_low = -BW/2;
    % f_high = BW/2;
    % b_lp = designfilt('lowpassfir', 'FilterOrder', 100, ...
    %     'CutoffFrequency', BW, 'SampleRate', Fs);
    % 
    % x_bb =  filtfilt(b_lp.Coefficients, 1,mixed_x);
    % freqz(b_lp.Coefficients,1,1024,Fs);
    
    % 2. Account for Group Delay
    % The total delay is span symbols (half from TX, half from RX if matched)
    
    [~, x_decoded] = pulseSHP(x_pb_noisy, simParams,'demodulate');
    message_rec = DiZeT(R, theta_c, x_decoded, K);
    
%     if figFlag
%         figure;
%         plot(1:9, [abs(x_decoded), abs(x)]);
%     end

    decoded_alphas = roots(flip(x_decoded'));
    %disp(decoded_alphas);
    
    %display decoded zeros
    if figFlag
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