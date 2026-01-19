% EbNo script

clear all;
close all;

% -- Params --
simParams.B = 1; %messages num
simParams.K = 400; %bits per message
simParams.Tsym = 1e-3;
simParams.fc = 15e3;
simParams.Fs = 200e3;
simParams.lambda = 1; %according to articles
simParams.figFlag = 0;



BW = 1 / simParams.Tsym; %assuming beta = 1
BWn = 3 * BW;
EbNo_dB = 0:18;
SNRa = EbNo_dB + 10*log10(BW / BWn) - 20; %in BMOCZ, 1 symbol is 1 bit
BER = zeros(length(EbNo_dB), 1);

for i = 1 : length(EbNo_dB)
    SNR = SNRa(i);
    [M_enc, M_dec] = MOCZstart(simParams, SNR);
    errosNum = sum(abs(M_enc - M_dec), 'all');
    BER(i) = errosNum / (simParams.B * simParams.K);
end

figure('Color', 'w'); % White background for reports
% Use semilogy for the logarithmic Y-axis
semilogy(EbNo_dB, BER, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
grid on;
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off'); % Improves readability
xlabel('E_b/N_0 (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('System BER Performance', 'FontSize', 14);
% Optional: Set specific limits to make it look "standard"
ylim([1e-6 1]); 
xlim([min(EbNo_dB) max(EbNo_dB)]);

legend('Simulated System', 'Location', 'southwest');

