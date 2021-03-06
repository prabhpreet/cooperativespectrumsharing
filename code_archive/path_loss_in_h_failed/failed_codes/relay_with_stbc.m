clear;

%%%%%%%%%
%PT Setup

pt_bits = 6;
pt_M = 2.^pt_bits;


%%%%%%%%%%%
%STBC Setup

rand('state',sum(100*clock));
randn('state',sum(100*clock));

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


%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set SNR here onwards for loop, when calculating BER
snr = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Notations

% = HX + N

%Y: M_r * T
%H: M_r * M_t
%X: M_t * T
%N: M_r * T

%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
%Phase I
%%%%%%%%

pt_symbols = randi([0 pt_M-1], 1,1)

pt_symbols = qammod(pt_symbols, pt_M, 0);

%At primary reciever
	pt_pr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

	pt_pr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr))*((randn(1,1)+i*randn(1,1))./sqrt(2));

	pt_pr_y = pt_pr_h*pt_symbols + pt_pr_n;

%At secondary transmitter
	
	pt_st_h = (randn(3,1)+i*randn(3,1))./sqrt(2);

	pt_st_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr))*((randn(3,1)+i*randn(3,1))./sqrt(2));

	pt_st_y = pt_st_h*pt_symbols + pt_st_n;

	pt_w = pt_st_h/norm(pt_st_h, 'fro'); %MRC beamformer
	
	pt_y_beamformed = pt_w' * pt_st_y;
	
	pt_st_y_decoded = qamdemod(pt_y_beamformed, pt_M, 0)
	
%At secondary receiver
	pt_sr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

	pt_sr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr))*((randn(1,1)+i*randn(1,1))./sqrt(2));

	pt_sr_y = pt_sr_h*pt_symbols + pt_sr_n;
	
	
%%%%%%%%
%Phase II
%%%%%%%%

%At secondary transmitter
	
	st_symbols = randi([1,4], 3, 1)
	st_symbols = st_su_const(st_symbols);

	c = st_pu_const(bi2de(reshape(de2bi(pt_st_y_decoded,pt_bits), 3, 2)) + ones(3,1));
	
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
	
	N = sqrt((norm(st_stbc_code,'fro')^2)/(4*snr))*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
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
		
		pr_decoded_symbol = Decoded_Symb{ind};
		
		pr_y_decoded = pr_decoded_symbol(1:3);
		
		pr_recieved=bi2de(reshape([(sign(imag(pr_y_decoded))+1)/2;(sign(real(pr_y_decoded))+1)/2].', 1,6))
		
		
	
%At secondary receiver
	h = (randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
		
	%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
	H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
	
	%the following complex equivalent channel is for the Embedded Diversity Code
	%for 3 TX
	H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];
	
	N = sqrt((norm(st_stbc_code,'fro')^2)/(4*snr))*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
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
		
		%BIT MAPPING
		sr_recieved = bi2de([(sign(imag(sr_y_decoded))+1)/2;(sign(real(sr_y_decoded))+1)/2].')+ones(3,1)