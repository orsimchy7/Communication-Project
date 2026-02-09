function [M, message_rec] = MOCZmodem(simParams)

    % -- Params --
    B = simParams.B; %num of sequences in a message
    L = simParams.L; % additional taps that the channel conv adds
    K = simParams.K ; %X(z) polynomial order. there are K+1 coefficents
    Tsym = simParams.Tsym;
    fc = simParams.fc;
    Fs = simParams.Fs;
    lambda = simParams.lambda;
    figFlag = simParams.figFlag;

    Lt = L + K; %number of coefficients in the recieved sequence. result of convolution with channel.
    zeroPad = zeros(L, B); %guard zeros - we allow the multipath components to exit the system.

    Kidxs = 1 : (K);
    
    M = randi([0,1],K ,B); %Binary Message.
    
    % -- Huffman BMOCZ --
    R = sqrt(1+2*lambda*sin(pi/K));
    theta_c = ((2*pi) * (Kidxs / K))';
    
    alphas = ((1-M)*(R^(-1)) + M * R).*exp(1i*theta_c);
    
    x_un = flip(poly(alphas))'; %finding coefficients and flip - now the free coefficent is first
    x = (x_un / norm(x_un)) * sqrt(K+1); % normalize to sqrt(K+1). now total energy of signal will be K+1
    x = [x ; zeroPad]; %Guard of zeros (channel of L taps)
    
    [x_pb, ~] = pulseSHP(x, simParams, 'modulate');

%     % -- Visualization --
    % if figFlag
    %     figure(1);
    %     subplot(2,1,1);
    %     stem(h_coeff, 'filled'); title('Raised Cosine Impulse Response'); grid on;
    %     subplot(2,1,2);
    %     t = 1:sps*(K+1); % Time axis for 9 symbols
    %     sig_segment = x_mod(t);
    %     % plot3(X, Y, Z) where X is time, Y is Real, Z is Imaginary
    %     plot3(t, real(sig_segment), imag(sig_segment), 'b', 'LineWidth', 1);
    %     hold on;
    %     % Add the discrete symbols as points in 3D space
    %     t_symbols = 1:sps:sps*(K+1);
    %     stem3(t_symbols, real(x), imag(x), 'r', 'LineWidth', 2);
    %     view(-25, 30); % Rotate for a better 3D perspective
    %     xlabel('Time (Samples)'); ylabel('Real (In-phase)'); zlabel('Imaginary (Quadrature)');
    %     title('3D Complex Pulse Shaped Signal');
    %     grid on;
    %     hold off;
    % end

    
    %% Decoding
    
    [~, x_decoded] = pulseSHP(x_pb, simParams,'demodulate');
    message_rec = DiZeT(R, theta_c, x_decoded, K);
    
%     if figFlag
%         figure;
%         plot(1:9, [abs(x_decoded), abs(x)]);
%     end

    decoded_alphas = roots(flip(x_decoded'));
    
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