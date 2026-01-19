function [signal_pb, x_decoded] = pulseSHP(x, simParams, mode)
%     x is:
%         1. in 'modulate': the payload of the packet
%         2. in 'demodulate': 
    
    % -- Params --
    L = simParams.L;
    K = simParams.K;
    N = K + L;
    Tsym = simParams.Tsym;
    fc = simParams.fc;
    Fs = simParams.Fs;
 
    sps = Tsym * Fs; %samples per symbol
    span = 6;
    beta = 1; % how much is beta?
    h_coeff = rcosdesign(beta, span, sps); %taps num is span*sps + 1
    normalization_factor = 1;
    group_delay = span * sps;
    t = (0:1/Fs:Tsym*(N+1 +group_delay/sps) - 1/Fs)';
    h_channel = [1; ...
             0.1982 - 0.975i; ...
             -0.809 - 0.421i; ...
             0.2 + 0.1i; ... % Placeholder for tap 4
             0.1 - 0.05i];   % Placeholder for tap 5


    if strcmp(mode, 'modulate')
        % -- Pulse Shaping --
        x_cor = conv(x, h_channel);
        x_upsamp = upsample(x_cor, sps);
        x_upsamp = [x_upsamp; zeros(group_delay,1)];
        signal_bb = filter(h_coeff, normalization_factor, x_upsamp);

%         h_channel_up = upsample(h_channel, sps);
%         signal_bb_faded = conv(h_channel_up, signal_bb);
        x_decoded = [];
        signal_pb = sqrt(2) * real(signal_bb.*exp(1i*2*pi*fc*t));
    elseif strcmp(mode, 'demodulate')
        % -- Decoding raised cosing --
        t = (0:1/Fs:Tsym*(K + 1 + (group_delay/sps)) - 1/Fs)';
        mixed_x = x.*exp(-1i*2*pi*fc*t);
        x_bb_filtered = filter(h_coeff, normalization_factor, mixed_x);
        
        % 3. Downsampling (The actual "Decoding" step)
        % We start sampling after the delay and then take 1 sample every 'sps'
        x_decoded = sqrt(2) * x_bb_filtered(group_delay + 1 : sps : end);
        signal_pb = [];
    end
end