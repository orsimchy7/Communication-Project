
% analizing results of generated MOCZ signal - BER
clear all;
close all;
%% load variables from mat file whose created by GenerateMOCZsignal
% load("L5B200K12N63u2.mat"); % 'signal_pb_total', 'P', 'simParams'
[~, ~, simParams, ~] = generateFunc();

%% preparing
B = simParams.B;
K = simParams.K;
BW = 1 / simParams.Tsym; %assuming beta = 1
BWn = 3 * BW;
EbNo_dB = 0:18;
% SNRa = EbNo_dB + 10*log10(BW / BWn); %in BMOCZ, 1 symbol is 1 bit
% SNRa = -10 : 1 : 6;
SNRa = 10*log10((K/((K+1)*simParams.Tsym*BW*(1+simParams.betha))) * 10.^(EbNo_dB/10));
% SNRa = ones(19,1)*100;
BER = zeros(length(SNRa), 1);
errorsNum = zeros(length(SNRa), 1);

attemptNum = 200 + 30 * (1:1:length(SNRa)); %200

%% running for different additive noise SNR
for i = 1 : length(SNRa)
    fprintf('i = %d \n', i);
    SNR = SNRa(i);
    for j = 1: attemptNum(i)
        [P, signal_pb_total, ~, x_bb, group_delay] = generateFunc();
        %[P_rec] = MOCZsimChannelNdecoding(simParams, SNR, signal_pb_total, P, x_bb, group_delay);

        %New Decoder (includes zadoffchu crosscorelation)
        [P_rec] = MOCZRealDecoding(simParams, SNR, signal_pb_total, P);
        errorsNum(i) = errorsNum(i) + sum(abs(P - P_rec), 'all');
    end
    BER(i) = errorsNum(i) / (B * K * attemptNum(i));
    fprintf('BER of snr %d is %d \n', SNR, BER(i));
end

%% plotting BER results
disp(errorsNum);
figure('Color', 'w'); % White background for reports
% Use semilogy for the logarithmic Y-axis
semilogy(EbNo_dB, BER, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
grid on;
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off'); % Improves readability
xlabel('E_b/N_0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('System BER Performance - normal channel and with zc', 'FontSize', 14);
% Optional: Set specific limits to make it look "standard"
% ylim([1e-6 1]); 
xlim([min(EbNo_dB) max(EbNo_dB)]);
legend('Simulated System', 'Location', 'southwest');
