%%  Plot graph
% Creates plot

clear;
close all;

snr_dB = -5:5:20;

files = ls('values*');
b = [];
x = [];
s = [];
for i = 1:size(files,1)
    load(files(i,:),'ber*');
    b = [b,ber_pu_direct];
    x = [x,ber_pu_relay];
    s = [s,ber_su];

end
ber_pu_direct = b;
ber_pu_relay = x;
ber_su = s;

semilogy(snr_dB, ber_pu_direct, 'r')
hold on
semilogy(snr_dB, ber_pu_relay, 'g')
hold on
semilogy(snr_dB, ber_su)
legend('BPSK Primary System', 'Relayed Primary System', 'Secondary System');
grid on;