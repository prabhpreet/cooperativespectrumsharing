%%  Compare graph
% Compares plots

clear all;
close all;

snr_dB = -5:5:30;
%snr_dB = snr_dB + 5;

[filename, pathname] = uigetfile('*.mat', 'Select data files', 'Multiselect', 'on')
files = strcat(char(pathname), filename);
b = [];
x = [];
s = [];
for i = 1:size(files, 2)
    load(files{i},'ber*', 'pt_bits');
    b = [b,ber_pu_direct];
    x = [x,ber_pu_relay];
    s = [s,ber_su];

end
ber_pu_direct_1 = b;
ber_pu_relay_1 = x;
ber_su_1 = s;

throughput_pu_direct_1 = pt_bits*(1-ber_pu_direct_1);
throughput_pu_relay_1 = pt_bits*(1-ber_pu_relay_1);
throughput_su_1 = 6*(1-ber_su_1);

[filename, pathname] = uigetfile('*.mat', 'Select data files', 'Multiselect', 'on')
files = strcat(char(pathname), filename);
b = [];
x = [];
s = [];
for i = 1:size(files, 2)
    load(files{i},'ber*', 'pt_bits');
    b = [b,ber_pu_direct];
    x = [x,ber_pu_relay];
    s = [s,ber_su];

end
ber_pu_direct_2 = b;
ber_pu_relay_2 = x;
ber_su_2 = s;

throughput_pu_direct_2 = pt_bits*(1-ber_pu_direct_2);
throughput_pu_relay_2 = pt_bits*(1-ber_pu_relay_2);
throughput_su_2 = 6*(1-ber_su_2);

figure
semilogy(snr_dB, ber_pu_direct_1, 'r')
hold on
semilogy(snr_dB, ber_pu_direct_2, 'g')
grid on;
title('BER PU DIRECT');
figure
semilogy(snr_dB, ber_pu_relay_1, 'r')
hold on
semilogy(snr_dB, ber_pu_relay_2, 'g')
grid on;
title('BER PU Relay');
figure
semilogy(snr_dB, ber_su_1, 'r')
hold on
semilogy(snr_dB, ber_su_2, 'g')
grid on;
title('BER SU');

figure
plot(snr_dB, throughput_pu_direct_1, 'r')
hold on
plot(snr_dB, throughput_pu_direct_2, 'g')
grid on
title('THROUGHPUT PU DIRECT');
figure
plot(snr_dB, throughput_pu_relay_1, 'r')
hold on
plot(snr_dB, throughput_pu_relay_1, 'g')
grid on
title('THROUGHPUT PU RELAY');
figure
plot(snr_dB, throughput_su_1,'r')
hold on
plot(snr_dB, throughput_su_1,'g')
grid on;
title('THROUGHPUT SU');