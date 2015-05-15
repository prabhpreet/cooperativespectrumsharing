%this code is for 3 TX and two layers of diversity, namely 3 and 1
%with rate 6/4


clear 

rand('state',sum(100*clock));
randn('state',sum(100*clock));

%QPSK symbols
Con_PSK4 = (exp(j.*[-3*pi/4 3*pi/4 7*pi/4 -7*pi/4])).*sqrt(2); %Unity amplitude and thus power

%scaled QPSK symbols
Con_PSK_red_pow = Con_PSK4./1.6;  %Scaling down the unity amplitude by 1.6

%create all combinations of two QPSK symbols: Why create these combinations?
a=[1:4 1:4 1:4];
b=unique([nchoosek(a,3)],'rows'); 
	%nchoosek: produces all combinations of vector where 3 elements are chosen at a time. unique: Choose unique combinations, since each element was assumed to be unique in combination.
	%Why choose 3 if creating combinations of two QPSK symbols
X = Con_PSK_red_pow(b((1:length(b)),:));

%just a naming convention
set1 = zeros(1,4);
set1 = Con_PSK4; %Message set A

set2 = zeros(1,4);
set2 = Con_PSK_red_pow; %Message set B

SNR = [5:5:25];
snr_lin=10.^(SNR/10);

err_count = 200; %[200 200 200 200 150];
snr_cnt = 1;

while snr_cnt <= length(snr_lin)
    
    snr = snr_lin(snr_cnt)
    Tot_Err_A = 0; %diversity 3 layer
    Tot_Err_B = 0; %diversity 1 layer
    no_packet(snr_cnt) = 1; %number of packets
  
	while Tot_Err_A < 200
		
		%RANDOM SYMBOLS
		x = randi([1 4],1,3); %3 QPSK symbols for diversity 3 layer
		c = set1(x);
		y = randi([1 4],1,3); %3 QPSK symbols for diversity 1 layer
		s = set2(y);
		
		Code = [c(1) c(2) c(3); -c(2)' c(1)' s(1)'; -c(3)' s(2)' c(1)'; s(3)' -c(3)' c(2)']; 
			%Conjugate: '
			%Row vectors: Instance of time, T= 4, Column vector: Tx Antennas M_t = 3;

		h = (randn(3,1)+i*randn(3,1))./sqrt(2); %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
		
		%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
		H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];
		
		%the following complex equivalent channel is for the Embedded Diversity Code
		%for 3 TX
		H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];

		%we use the average energy per symbol duration in the codeword
		N = sqrt((norm(Code,'fro')^2)/(4*snr))*((randn(4,1)+i*randn(4,1))./sqrt(2)); 
			% T= 4, M_r = 1
			%Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.
		
		%received Signal
		Y = Code * h  + N; %Code multiplied with h, not H i.e. multiplied with fading coefficients only. 
		
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
		
		MLDecod_Symb = Decoded_Symb{ind};
		
		%BIT MAPPING
		Rx_bit = [];
		Rx_bit(:,1:2:12)=(sign(real(MLDecod_Symb))+1)/2;
		Rx_bit(:,2:2:12)=(sign(imag(MLDecod_Symb))+1)/2;
		Rx_bit = transpose(Rx_bit);
		
		Tx_bit=[];
		Tx_bit(:,1:2:12)=(sign(real([c s]))+1)/2;
		Tx_bit(:,2:2:12)=(sign(imag([c s]))+1)/2;
		Tx_bit = transpose(Tx_bit);

		
		%BIT ERROR RATE
		
		Error_A = (Rx_bit(1:6,:)~=Tx_bit(1:6,:))';
		Error_B = (Rx_bit(7:12,:)~=Tx_bit(7:12,:))';
	   
		Tot_Err_A = Tot_Err_A + sum(Error_A);
		Tot_Err_B = Tot_Err_B + sum(Error_B);
		
		if mod(no_packet(snr_cnt),100000)==0 %Display total errors for every 10,000 packets.
			no_packet(snr_cnt)
			Tot_Err_A
		end
		no_packet(snr_cnt) = no_packet(snr_cnt)+1;
	end
       
    BER_A(snr_cnt) =  Tot_Err_A/(6*no_packet(snr_cnt)) %3 symbols are sent per packet per message set. Therefore, 6 bits per packet.
    BER_B(snr_cnt) =  Tot_Err_B/(6*no_packet(snr_cnt))
    
    snr_cnt = snr_cnt + 1;
    
end

%save NewCode_4PSK_Rate6by4_HybMLIC_AvgSNR_K2 BER_A BER_B SNR
%Plot the SNR vs BER curve

semilogy(SNR,BER_A,'bo-', SNR,BER_B,'r*-')
    
xlabel('SNR in dB');
ylabel('Bit Error Rate');
axis([5,30,10^(-7),10^0]);
grid
%legend('BER for Diversity Layer 2','BER for Diversity Layer 3');
    
   
    

 