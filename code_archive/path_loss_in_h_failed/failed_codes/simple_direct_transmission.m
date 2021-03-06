clear;
close;

rng('default');
rng('shuffle');

snr_dB = 5:5:30
snr = 10.^(snr_dB./10);

bpsk_constellation = [-1 1];

ber_pu_direct = zeros(1,length(snr));

for i = 1:length(snr)

	errors_pu_direct = 0;
	no_packet_pu_direct = 0;
	while errors_pu_direct < 3000
		pt_direct_bits = randi([0 1], 1, 6);
		pt_direct_x = bpsk_constellation(pt_direct_bits + 1);
		
		pt_direct_h = randn(1, 1) + j.*randn(1,1)./sqrt(2);
		
		pt_direct_n = sqrt((norm(pt_direct_x(1),'fro')^2)/(snr(i)))*((randn(1,1)+j*randn(1,1))./sqrt(2));
		
		pt_direct_y = pt_direct_h*pt_direct_x  + pt_direct_n;
		
		%Perfect CSI
		
		pt_direct_y = pt_direct_h'*pt_direct_y./(norm(pt_direct_h, 'fro').^2);
		
		pt_direct_decoded = (sign(real(pt_direct_y))+1)./2;
		
		errors_pu_direct = errors_pu_direct + sum(pt_direct_bits ~= pt_direct_decoded);

		no_packet_pu_direct = no_packet_pu_direct + 1;
	end
	
	ber_pu_direct(i) = errors_pu_direct./no_packet_pu_direct;
end


grid on;
semilogy(snr_dB, ber_pu_direct);
legend('Direct BSPK from Primary Tx to Rx');