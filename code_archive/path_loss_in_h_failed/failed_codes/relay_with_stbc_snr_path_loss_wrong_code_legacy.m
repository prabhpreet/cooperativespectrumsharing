clear;
close;

%%%%%%%%%
%System model

d_pt_pr = 1;
d_pt_st = 0.5;
d_st_pr = 0.5;
d_st_sr = 0.25;

m = -3;

d_pt_pr = d_pt_pr.^m;
d_pt_st = d_pt_st.^m;
d_st_pr = d_st_pr.^m;
d_st_sr = d_st_sr.^m;

%%%%%%%%%
%PT Setup

pt_bits = 6;
pt_M = 2.^pt_bits;



%%%%%%%%%%%
%STBC Setup

rng('default');
rng('shuffle');

%QPSK symbols
st_pu_const = (exp(j.*[-3*pi/4 3*pi/4 7*pi/4 -7*pi/4])).*sqrt(2); %Unity magnitude and thus unity power

%scaled QPSK symbols
st_su_const = st_pu_const./1.6;  %Scaling down the unity amplitude by 1.6

%create all combinations of two QPSK symbols
a=[1:4 1:4 1:4];
b=unique([nchoosek(a,3)],'rows'); 
	%nchoosek: produces all combinations of vector where 3 elements are chosen at a time. unique: Choose unique combinations, since each element was assumed to be unique in combination.
X = st_su_const(b((1:length(b)),:));


%%%%%%%%%%%

snr_dB = 5:5:30;

snr = 10.^(snr_dB./10);

ber_pu_direct = zeros(1,length(snr));
ber_pu_relay = zeros(1,length(snr));
ber_su = zeros(1,length(snr));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set SNR here onwards for loop, when calculating BER
for m = 1:length(snr)

	no_tx_pu = 0;
	no_tx_su = 0;
	errors_pu_direct = 0;
	errors_pu_relay = 0;
	errors_su = 0;

	while errors_pu_direct < 200
	

		%%%%%%%%%%%%%%%%%%%%%%%%%
		%Notations

		%Y = HX + N (except for STBC, which is Y=XH + N)

		%Y: M_r * T
		%H: M_r * M_t
		%X: M_t * T
		%N: M_r * T

		%%%%%%%%%%%%%%%%%%%%%%%%%

		%%%%%%%%
		%Direct Transmission of PU
		%%%%%%%%
		
		pt_direct_symbols = randi([0 pt_M-1], 1,1);
		
		pt_direct_bits = de2bi(pt_direct_symbols, pt_bits);

		pt_direct_bit_sequence = de2bi(pt_direct_symbols,pt_bits);

		pt_direct_symbols = qammod(pt_direct_symbols, pt_M, 0);

		%At primary reciever
			pt_direct_pr_h = d_pt_pr.*(randn(1,1)+i*randn(1,1))./sqrt(2);

			pt_direct_pr_n = sqrt((norm(pt_direct_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(1,1)+i*randn(1,1))./sqrt(2));

			pt_direct_pr_y = pt_direct_pr_h*pt_direct_symbols + pt_direct_pr_n;
			
			%PERFECT CSI
			pt_direct_pr_y = pt_direct_pr_h'*pt_direct_pr_y./(norm(pt_direct_pr_h, 'fro').^2);
			
			pt_direct_pr_y_decoded = qamdemod(pt_direct_pr_y, pt_M, 0);
			
			pt_direct_pr_y_decoded_bits = de2bi(pt_direct_pr_y_decoded,pt_bits);
			
			errors_pu_direct = errors_pu_direct + sum(pt_direct_bits ~= pt_direct_pr_y_decoded_bits);
		
		
		%%%%%%%%
		%Phase I
		%%%%%%%%

		pt_symbols = randi([0 pt_M-1], 1,1);

		pt_bit_sequence = de2bi(pt_symbols,pt_bits);

		pt_symbols = qammod(pt_symbols, pt_M, 0);

		%At primary reciever
			pt_pr_h =  d_pt_pr.*(randn(1,1)+i*randn(1,1))./sqrt(2);

			pt_pr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(1,1)+i*randn(1,1))./sqrt(2));

			pt_pr_y = pt_pr_h*pt_symbols + pt_pr_n;
			
			pt_pr_y = pt_pr_h'*pt_pr_y ./(norm(pt_pr_h, 'fro').^2);

		%At secondary transmitter
			
			pt_st_h = d_pt_st.*(randn(3,1)+i*randn(3,1))./sqrt(2);

			pt_st_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(3,1)+i*randn(3,1))./sqrt(2));

			pt_st_y = pt_st_h*pt_symbols + pt_st_n;
			
			%%%%%%%
			%Doubt: Beamforming and channel estimation? Divide channel estimate before beamforming? Is effect of increasing recieve diversity significant after channel estimation techniques?
			
			pt_w = pt_st_h/norm(pt_st_h, 'fro'); %MRC beamformer
			
			%pt_y_beamformed = pt_w' * pt_st_y; %DOUBT!
			
			%%%%%%
			%Temporary: Assume one recieve antenna only!
			pt_y_beamformed = pt_st_h(1)'.*pt_st_y(1)./(norm(pt_st_h(1), 'fro').^2);
			%%%%%%
				
			pt_st_y_decoded = qamdemod(pt_y_beamformed, pt_M, 0);
			
			pt_st_y_decoded_bits = de2bi(pt_st_y_decoded,pt_bits);
			
		%At secondary receiver
			pt_sr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

			pt_sr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(1,1)+i*randn(1,1))./sqrt(2));

			pt_sr_y = pt_sr_h*pt_symbols + pt_sr_n;
			
			
			
		%%%%%%%%
		%Phase II
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
			h = d_st_pr.*(randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
			%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
			H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
			
			%the following complex equivalent channel is for the Embedded Diversity Code
			%for 3 TX
			H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];
			
			N = sqrt((norm(st_stbc_code,'fro')^2)/(4*snr(m)))*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
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
			h = d_st_sr.*(randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
				
			%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
			H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
			
			%the following complex equivalent channel is for the Embedded Diversity Code
			%for 3 TX
			H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];
			
			N = sqrt((norm(st_stbc_code,'fro')^2)/(4*snr(m)))*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
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
		
		snr(m)
		ber_pu_direct(m) = errors_pu_direct/(pt_bits*no_tx_pu)
		ber_pu_relay(m) = errors_pu_relay/(pt_bits*no_tx_pu)
		ber_su(m) = errors_su/(6*no_tx_su)
end

semilogy(snr_dB, ber_pu_direct, 'r')
hold on
semilogy(snr_dB, ber_pu_relay, 'g')
hold on
semilogy(snr_dB, ber_su)
legend('64 QAM Direct Primary System', 'Relayed Primary System', 'Secondary System');
grid on;