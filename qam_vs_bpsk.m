clear;
close;

%%%%%%%%%
%Pathloss model: d^-n, d is normalized distance

% Let's assume all are collinear. PT---ST---SR----PR

d_pt_pr = 10; %in meters!
d_pt_st = 0.8;
d_st_pr = d_pt_pr - d_pt_st;
d_st_sr = 5;

n = 2; %Free space path loss



%%%%%%%%%
%Primary Transmitter Setup

pt_direct_constellation = [-1 1];
pt_direct_average_symbol_energy = mean(abs(pt_direct_constellation).^2);

%Relay Phase
pt_bits = 6;
pt_M = 2.^pt_bits;
pt_relay_constellation = qammod([0 pt_M-1], pt_M, 0);
pt_relay_average_symbol_energy = mean(abs(pt_relay_constellation).^2);


%%%%%%%%%%%
%STBC Setup

rng('default');
rng('shuffle');

%QPSK symbols
st_pu_const = (exp(j.*[-3*pi/4 3*pi/4 7*pi/4 -7*pi/4])).*sqrt(2); %Unity magnitude and thus unity power
st_pu_average_symbol_energy = mean(abs(st_pu_const).^2);

%scaled QPSK symbols
st_su_const = st_pu_const./1.6;  %Scaling down the unity amplitude by 1.6
st_su_average_symbol_energy = mean(abs(st_su_const).^2);

st_timeslot_energy = [3*st_pu_average_symbol_energy; (2*st_pu_average_symbol_energy) + st_su_average_symbol_energy; (2*st_pu_average_symbol_energy) + st_su_average_symbol_energy;(2*st_pu_average_symbol_energy) + st_su_average_symbol_energy];


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
	
	pt_direct_pr_snr = snr_dB(m) - (n.*10.*log10(d_pt_pr));
	pt_direct_pr_snr = 10.^(pt_direct_pr_snr/10);
	pt_direct_pr_sigma = sqrt(pt_direct_average_symbol_energy/(4.*pt_direct_pr_snr));
	
	pt_relay_pr_snr = snr_dB(m) - (n.*10.*log10(d_pt_pr));
	pt_relay_pr_snr = 10.^(pt_relay_pr_snr/10);
	pt_relay_pr_sigma = sqrt(pt_relay_average_symbol_energy/(4.*pt_relay_pr_snr));
	
	pt_st_snr = snr_dB(m) - (n.*10.*log10(d_pt_st));
	pt_st_snr = 10.^(pt_st_snr/10);
	pt_st_sigma = sqrt(pt_relay_average_symbol_energy/(4.*pt_st_snr));
	
	
	st_pr_snr = snr_dB(m) - (n.*10.*log10(d_st_pr));
	st_pr_snr = 10.^(st_pr_snr/10);
	st_pr_sigma = sqrt(st_timeslot_energy./(4.*st_pr_snr));
	%st_pr_sigma = (sqrt([norm(st_stbc_code(1,:),'fro') ;norm(st_stbc_code(2,:),'fro');norm(st_stbc_code(3,:),'fro');norm(st_stbc_code(4,:),'fro')]).^2)./(2.*st_pr_snr);
	
	
	st_sr_snr = snr_dB(m) - (n.*10.*log10(d_st_sr));
	st_sr_snr = 10.^(st_sr_snr/10);
	st_sr_sigma = sqrt(st_timeslot_energy./(4.*st_sr_snr));
	%st_sr_sigma = (sqrt([norm(st_stbc_code(1,:),'fro') ;norm(st_stbc_code(2,:),'fro');norm(st_stbc_code(3,:),'fro');norm(st_stbc_code(4,:),'fro')]).^2)./(2.*st_sr_snr);

	while errors_pu_direct < 1000
	

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
		%BSPK, 6 bits at a time. h constant for whole time interval.
		%%%%%%%%
		
		
		pt_direct_bits = randi([0 1], 1, 6);
		pt_direct_x = pt_direct_constellation(pt_direct_bits + 1);
		
		pt_direct_h = randn(1, 1) + j.*randn(1,1)./sqrt(2);
		
		pt_direct_n = pt_direct_pr_sigma*((randn(1,1)+i*randn(1,1))./sqrt(2));
		
		pt_direct_y = pt_direct_h*pt_direct_x  + pt_direct_n;
		
		%Perfect CSI
		
		pt_direct_y = pt_direct_h'*pt_direct_y./(norm(pt_direct_h, 'fro').^2);
		
		pt_direct_decoded = (sign(real(pt_direct_y))+1)./2;
		
		errors_pu_direct = errors_pu_direct + sum(pt_direct_bits ~= pt_direct_decoded);
		
		
		%%%%%%%%
		%Relayed System
		%Phase I: PT, 64QAM since PT near ST (relay).
		%%%%%%%%

		pt_symbols = randi([0 pt_M-1], 1,1);

		pt_bit_sequence = de2bi(pt_symbols,pt_bits);

		pt_symbols = qammod(pt_symbols, pt_M, 0);

		%At primary reciever
			pt_pr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

			pt_pr_n = pt_relay_pr_sigma*((randn(1,1)+i*randn(1,1))./sqrt(2));

			pt_pr_y = pt_pr_h*pt_symbols + pt_pr_n;
			
			pt_pr_y = pt_pr_h'*pt_pr_y ./(norm(pt_pr_h, 'fro').^2);

		%At secondary transmitter
			
			pt_st_h = (randn(3,1)+i*randn(3,1))./sqrt(2);

			pt_st_n = pt_st_sigma*((randn(3,1)+i*randn(3,1))./sqrt(2));

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
			
			errors_pu_relay = errors_pu_relay + sum(pt_bit_sequence ~= pt_st_y_decoded_bits);
		

			
		% %At secondary receiver
			% pt_sr_h = (randn(1,1)+i*randn(1,1))./sqrt(2);

			% pt_sr_n = sqrt((norm(pt_symbols,'fro')^2)/(pt_bits*snr(m)))*((randn(1,1)+i*randn(1,1))./sqrt(2));

			% pt_sr_y = pt_sr_h*pt_symbols + pt_sr_n;
			
			no_tx_pu = no_tx_pu + 1;
			no_tx_su = no_tx_su + 1;
		
		end
		
		snr_dB(m)
		ber_pu_direct(m) = errors_pu_direct/(pt_bits*no_tx_pu)
		ber_pu_relay(m) = errors_pu_relay/(pt_bits*no_tx_pu)
		ber_su(m) = errors_su/(6*no_tx_su)
end

semilogy(snr_dB, ber_pu_direct, 'r')
hold on
semilogy(snr_dB, ber_pu_relay, 'g')
hold on
semilogy(snr_dB, ber_su)
legend('BPSK Primary System', 'Relayed Primary System', 'Secondary System');
grid on;