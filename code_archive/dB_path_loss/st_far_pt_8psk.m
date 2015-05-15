clear;
close;

%%%%%%%%%
%Pathloss model

% Let's assume all are collinear. PT---ST---SR----PR

pl_pt_pr_dB = 0;
pl_pt_st_dB = 5;
pl_st_pr_dB = 10;
pl_st_sr_dB = 0;

%%%%%%%%%
%Primary Transmitter Setup

pt_direct_constellation = [-1 1];
pt_direct_average_symbol_energy = mean(abs(pt_direct_constellation).^2);
pt_direct_constellation = pt_direct_constellation./sqrt(pt_direct_average_symbol_energy);
pt_direct_average_symbol_energy = mean(abs(pt_direct_constellation).^2);

%Relay Phase
pt_bits = 3;
pt_M = 2.^pt_bits;
pt_relay_constellation = pskmod([0:pt_M-1], pt_M);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);
pt_relay_constellation = pt_relay_constellation./sqrt(pt_relay_average_symbol_energy);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);



%%%%%%%%%%%
%STBC Setup

rng('default');
rng('shuffle');

%QPSK symbols
st_pu_const = (exp(j.*[-3*pi/4 3*pi/4 7*pi/4 -7*pi/4])); %Unity magnitude and thus unity power
st_pu_average_symbol_energy = mean(abs(st_pu_const).^2);


%scaled QPSK symbols
st_su_const = st_pu_const./1.6;  %Scaling down the unity amplitude by 1.6
st_su_average_symbol_energy = mean(abs(st_su_const).^2);

e = 4./((9*st_pu_average_symbol_energy) + (3*st_su_average_symbol_energy));
st_pu_const = st_pu_const.*sqrt(e);
st_su_const = st_su_const.*sqrt(e);



%create all combinations of two QPSK symbols
a=[1:4 1:4 1:4];
b=unique([nchoosek(a,3)],'rows'); 
	%nchoosek: produces all combinations of vector where 3 elements are chosen at a time. unique: Choose unique combinations, since each element was assumed to be unique in combination.
X = st_su_const(b((1:length(b)),:));


%%%%%%%%%%%

snr_dB = 5:5:30;


ber_pu_direct = zeros(1,length(snr_dB));
ber_pu_relay = zeros(1,length(snr_dB));
ber_su = zeros(1,length(snr_dB));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set SNR here onwards for loop, when calculating BER
for m = 1:length(snr_dB)

	no_tx_pu = 0;
	no_tx_su = 0;
	errors_pu_direct = 0;
	errors_pu_relay = 0;
	errors_su = 0;
	
	pt_direct_pr_snr = snr_dB(m) + pl_pt_pr_dB;
	pt_direct_pr_snr = 10.^(pt_direct_pr_snr/10);
	pt_direct_pr_sigma = sqrt(1/pt_direct_pr_snr);
	
	
	pt_st_snr = snr_dB(m)+ pl_pt_st_dB ;
	pt_st_snr = 10.^(pt_st_snr/10);
	pt_st_sigma = sqrt(1/pt_st_snr);
	
	
	st_pr_snr = snr_dB(m) +pl_st_pr_dB;
	st_pr_snr = 10.^(st_pr_snr/10);
	st_pr_sigma = sqrt(1/st_pr_snr);
	%st_pr_sigma = (sqrt([norm(st_stbc_code(1,:),'fro') ;norm(st_stbc_code(2,:),'fro');norm(st_stbc_code(3,:),'fro');norm(st_stbc_code(4,:),'fro')]).^2)./(2.*st_pr_snr);
	
	
	st_sr_snr = snr_dB(m) + pl_st_sr_dB;
	st_sr_snr = 10.^(st_sr_snr/10);
	st_sr_sigma = sqrt(1/st_sr_snr);
	%st_sr_sigma = (sqrt([norm(st_stbc_code(1,:),'fro') ;norm(st_stbc_code(2,:),'fro');norm(st_stbc_code(3,:),'fro');norm(st_stbc_code(4,:),'fro')]).^2)./(2.*st_sr_snr);

	while errors_pu_relay < 200
	

		%%%%%%%%%%%%%%%%%%%%%%%%%
		%Notations
		
		%Assumes h remains same for 6 T

		%Y = HX + N (except for STBC, which is Y=XH + N)

		%Y: M_r * T
		%H: M_r * M_t
		%X: M_t * T
		%N: M_r * T

		%%%%%%%%%%%%%%%%%%%%%%%%%

		%%%%%%%%
		%Direct Transmission of PU
		%BSPK, 6 bits at a time. h constant for whole time interval.
		%%%%%%%%
		
				
		pt_direct_bits = randi([0 1], 1, 6);
		pt_direct_x = pt_direct_constellation(pt_direct_bits + 1);
		
		pt_direct_h = randn(1, 6) + j.*randn(1,6)./sqrt(2);
		
		pt_direct_n = pt_direct_pr_sigma*((randn(1,6)+i*randn(1,6))./sqrt(2));
		
		pt_direct_y = pt_direct_h.*pt_direct_x  + pt_direct_n;
		
		%Perfect CSI
		
		pt_direct_y = conj(pt_direct_h).*pt_direct_y./(abs(pt_direct_h).^2);
		
		pt_direct_decoded = (sign(real(pt_direct_y))+1)./2;
		
		errors_pu_direct = errors_pu_direct + sum(pt_direct_bits ~= pt_direct_decoded);

		%%%%%%%%
		%Relayed System
		%Phase I: PT, 8PSK * 2 since PT near ST (relay).
		%%%%%%%%

		pt_symbols = randi([0 pt_M-1], 1,2);

		pt_bit_sequence = de2bi(pt_symbols,pt_bits);
		
		pt_bit_sequence = pt_bit_sequence(:)';

		pt_symbols = pskmod(pt_symbols, pt_M);


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
			
			pt_st_y_decoded_bits = pt_st_y_decoded_bits(:)';
		% %At secondary receiver
			% pt_sr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

			% pt_sr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(1,1)+i*randn(1,1))./sqrt(2));

			% pt_sr_y = pt_sr_h*pt_symbols + pt_sr_n;
			
			
		
		%%%%%%%%
		%Phase II: ST, STBC Code
		%%%%%%%%

		%At secondary transmitter
			
			
			
			st_symbols = randi([1,4], 3, 1);
			st_bit_sequence = de2bi(st_symbols-1, 2);
			st_symbols = st_su_const(st_symbols);

			
			c = st_pu_const(bi2de(reshape(pt_st_y_decoded_bits, 3, 2)) + ones(3,1));
			
			s = st_symbols;

			st_stbc_code = [c(1) c(2) c(3); -c(2)' c(1)' s(1)'; -c(3)' s(2)' c(1)'; s(3)' -c(3)' c(2)'];
					%Conjugate: '
					%Row vectors: Instance of time, T= 4, Column vector: Tx Antennas M_t = 3;
			
			
			
			
		%At primary receiver
			h = (randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
			%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
			H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
			
			%the following complex equivalent channel is for the Embedded Diversity Code
			%for 3 TX
			H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];
			
			N = st_pr_sigma.*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
					% T= 4, M_r = 1
					%Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
				%received Signal
				Y = st_stbc_code * h  + N; %Code multiplied with h, not H i.e. multiplied with fading coefficients only. 
				
				for k = 1:length(X)
					%discard the effect of diversity 2 and 1 layer from total received
					%signal to get Y_remaining, alias Y_rem
					Y_rem = [Y(1); Y(2)-h(3)*X(k,1)'; Y(3)-h(2)*X(k,2)'; Y(4)-h(1)*X(k,3)'];
					
					Y_prime = [Y_rem(1) Y_rem(2:4)'].';
					
					%Apply matched filtering because the remaining received signal is
					%due to the contribution from the Orthogonal Diversity 3 layer Only
					Y_match = H' * Y_prime;
					
					Sym = sign([real(Y_match); imag(Y_match)]);
				
					S_tilde = Sym(1:3) + i*Sym(4:6);
					
					Decoded_Symb{k} = [S_tilde.' X(k,1) X(k,2) X(k,3)];
					%Now apply ML decoding using the overall equivalent channel matrix
					%H_Eqv
					diff = [Y(1) Y(2:4)'].' - H_eqv * Decoded_Symb{k}.';
					
					metric(k) = norm(diff,'fro')^2;
					
				end
				
				[W, ind] = min(metric);
				
				st_pr_decoded_stbc = Decoded_Symb{ind};
				
				st_pr_decoded_symbols = st_pr_decoded_stbc(1:3);
				
				st_pr_decoded_bits = reshape([(sign(imag(st_pr_decoded_symbols))+1)/2;(sign(real(st_pr_decoded_symbols))+1)/2].', 1,6);
				
				pr_recieved=bi2de(st_pr_decoded_bits);
				
				errors_pu_relay = errors_pu_relay + sum(pt_bit_sequence ~= st_pr_decoded_bits);

				no_tx_pu = no_tx_pu + 1;
				
				
			
		%At secondary receiver
			h = (randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
			%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
			H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
			
			%the following complex equivalent channel is for the Embedded Diversity Code
			%for 3 TX
			H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];
			
			N = st_sr_sigma.*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
					% T= 4, M_r = 1
					%Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
				%received Signal
				Y = st_stbc_code * h  + N; %Code multiplied with h, not H i.e. multiplied with fading coefficients only. 
				
				for k = 1:length(X)
					%discard the effect of diversity 2 and 1 layer from total received
					%signal to get Y_remaining, alias Y_rem
					Y_rem = [Y(1); Y(2)-h(3)*X(k,1)'; Y(3)-h(2)*X(k,2)'; Y(4)-h(1)*X(k,3)'];
					
					Y_prime = [Y_rem(1) Y_rem(2:4)'].';
					
					%Apply matched filtering because the remaining received signal is
					%due to the contribution from the Orthogonal Diversity 3 layer Only
					Y_match = H' * Y_prime;
					
					Sym = sign([real(Y_match); imag(Y_match)]);
				
					S_tilde = Sym(1:3) + i*Sym(4:6);
					
					Decoded_Symb{k} = [S_tilde.' X(k,1) X(k,2) X(k,3)];
					%Now apply ML decoding using the overall equivalent channel matrix
					%H_Eqv
					diff = [Y(1) Y(2:4)'].' - H_eqv * Decoded_Symb{k}.';
					
					metric(k) = norm(diff,'fro')^2;
					
				end
				
				[W, ind] = min(metric);
				
				sr_decoded_symbol = Decoded_Symb{ind};

				sr_y_decoded = sr_decoded_symbol(4:6);
				
				sr_y_decoded_bits = [(sign(imag(sr_y_decoded))+1)/2;(sign(real(sr_y_decoded))+1)/2].';
				
				%BIT MAPPING
				sr_recieved = bi2de(sr_y_decoded_bits)+ones(3,1);
				
				no_tx_su = no_tx_su + 1;
				
				errors_su = errors_su + sum(st_bit_sequence(:)~=sr_y_decoded_bits(:));
		
		end
		
		snr_dB(m)
		ber_pu_direct(m) = errors_pu_direct/(pt_bits*no_tx_pu)
		ber_pu_relay(m) = errors_pu_relay/(pt_bits*no_tx_pu)
		ber_su(m) = errors_su/(6*no_tx_su)
		
end

figure
semilogy(snr_dB, ber_pu_direct, 'r')
hold on
semilogy(snr_dB, ber_pu_relay, 'g')
hold on
semilogy(snr_dB, ber_su)
legend('BPSK Primary System', 'Relayed Primary System', 'Secondary System');
grid on;