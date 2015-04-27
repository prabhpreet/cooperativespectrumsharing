%% Notations
%$$Y = HX + N$$
%(except for STBC, which is $Y=XH + N$)
%Y: M_r * T
%H: M_r * M_t
%X: M_t * T
%N: M_r * T

%% Initialization

clear;
close;

%Specify extra stings to save .mat file as (leave no spaces)
version = 'secondarySystemOnly';

%% SNR and error parameters
% Set SNR and errors to be evaulated

errors_evaluated = 100;

snr_dB = 5:5:30;

ber_su = zeros(1,length(snr_dB));

%% Pathloss Model
% Taking all pathlosses in dB

pl_pt_pr_dB = 0;
%pl_pt_st_dB = 10;
pl_st_pr_dB = 0;
pl_st_sr_dB = 0;

%% Primary Transmitter Setup

pt_direct_constellation = [-1 1];
pt_direct_average_symbol_energy = mean(abs(pt_direct_constellation).^2);
pt_direct_constellation = pt_direct_constellation./sqrt(pt_direct_average_symbol_energy);
pt_direct_average_symbol_energy = mean(abs(pt_direct_constellation).^2);


%Relay Phase: TODO- check
pt_bits = 6;
pt_M = 2.^pt_bits;
pt_relay_constellation = qammod([0:pt_M-1], pt_M, 0);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);
pt_relay_constellation = pt_relay_constellation./sqrt(pt_relay_average_symbol_energy);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);

%% STBC Setup

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


%% Loop section
% Set SNR here onwards for loop, when calculating BER

for m = 1:length(snr_dB)

	no_tx_su = 0;
	errors_su = 0;

	st_sr_snr = snr_dB(m) +pl_st_sr_dB;
	st_sr_snr = 10.^(st_sr_snr/10);
	st_sr_sigma = sqrt(1/st_sr_snr);
    
	while errors_su < errors_evaluated

		%%%%%%%%
		%Relayed System
		%Phase I: Assume perfect decoding.
		%%%%%%%%

		pt_symbols = randi([0 pt_M-1], 1,1);

		pt_bit_sequence = de2bi(pt_symbols,pt_bits);

		%At secondary transmitter
			
		pt_st_y_decoded = pt_symbols;
			
		pt_st_y_decoded_bits = de2bi(pt_st_y_decoded,pt_bits);
		
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
			
        % At secondary reciever
            st_sr_h = (randn(3,1)+1i*randn(3,1))./sqrt(2);
            st_sr_n = st_sr_sigma.*((randn(4,1)+1i*randn(4,1))./sqrt(2));
            
            st_sr_Y = st_stbc_code * st_sr_h+ st_sr_n;
            
            %%%%%%%%%%%%%%%
            % Assume that we know all 6 bits of primary transmitter
            % perfectly.
            %%%%%%%%%%%%%%%      
            
            st_sr_Y_dash = st_sr_Y(2:4);
            H_dash = [-st_sr_h(2) st_sr_h(1) 0;-st_sr_h(3) 0 st_sr_h(1); 0 -st_sr_h(3) st_sr_h(2)];
            H_prime = [st_sr_h(3)'./(abs(st_sr_h(3)).^2);st_sr_h(2)'./(abs(st_sr_h(2)).^2);st_sr_h(1)'./(abs(st_sr_h(1)).^2)];
            
            s_prime = (H_prime.*(st_sr_Y_dash + H_dash*(c')))';
            
            st_sr_decoded_bits = [(sign(imag(s_prime))+1)/2;(sign(real(s_prime))+1)/2].';
            
            
            errors_su = errors_su + sum(st_bit_sequence(:) ~= st_sr_decoded_bits(:))
            
            no_tx_su = no_tx_su + 1;
            
    end
		
		snr_dB(m)
		ber_su(m) = errors_su/(6*no_tx_su)
end

%% Save Data

save(strcat('values_', num2str(snr_dB(1)), '_' , num2str(snr_dB(end)), 'dB','_',num2str(errors_evaluated),'iterations', version,'.mat'));
