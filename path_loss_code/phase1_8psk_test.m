pt_bits = 3;
pt_M = 2.^pt_bits;
pt_relay_constellation = pskmod([0 pt_M-1], pt_M);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);

m = 1;
d_pt_pr = 1;
d_pt_st = 0.9;
n = 3;

snr_dB = 0.1:30;

	pt_relay_pr_snr = snr_dB(m) - (n.*10.*log10(d_pt_pr));
	pt_relay_pr_snr = 10.^(pt_relay_pr_snr/10);
	pt_relay_pr_sigma = sqrt(pt_relay_average_symbol_energy/(4.*pt_relay_pr_snr));
	
		pt_st_snr = snr_dB(m) - (n.*10.*log10(d_pt_st));
	pt_st_snr = 10.^(pt_st_snr/10);
	pt_st_sigma = sqrt(pt_relay_average_symbol_energy/(4.*pt_st_snr));

		pt_symbols = randi([0 pt_M-1], 1,2)

		pt_bit_sequence = de2bi(pt_symbols,pt_bits);
		
		pt_bit_sequence = pt_bit_sequence(:)'

		pt_symbols = pskmod(pt_symbols, pt_M);

		%At primary reciever
			pt_pr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

			pt_pr_n = pt_relay_pr_sigma*((randn(1,2)+i*randn(1,2))./sqrt(2));

			pt_pr_y = pt_pr_h*pt_symbols + pt_pr_n;
			
			pt_pr_y = pt_pr_h'*pt_pr_y ./(norm(pt_pr_h, 'fro').^2);

		%At secondary transmitter
			
			pt_st_h = (randn(3,1)+i*randn(3,1))./sqrt(2);

			pt_st_n = pt_st_sigma*((randn(3,2)+i*randn(3,2))./sqrt(2));

			pt_st_y = pt_st_h*pt_symbols + pt_st_n;
			
			%%%%%%%
			%Doubt: Beamforming and channel estimation? Divide channel estimate before beamforming? Is effect of increasing recieve diversity significant after channel estimation techniques?
			
			pt_w = pt_st_h/norm(pt_st_h, 'fro'); %MRC beamformer
			
			%pt_y_beamformed = pt_w' * pt_st_y; %DOUBT!
			
			%%%%%%
			%Temporary: Assume one recieve antenna only!
			pt_y_beamformed = pt_st_h(1,:)'.*pt_st_y(1,:)./(norm(pt_st_h(1,:), 'fro').^2);
			%%%%%%
				
			pt_st_y_decoded = pskdemod(pt_y_beamformed, pt_M);
			
			pt_st_y_decoded_bits = de2bi(pt_st_y_decoded,pt_bits);
			
			pt_st_y_decoded_bits = pt_st_y_decoded_bits(:)'