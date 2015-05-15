pt_bits = 6;
pt_M = 2.^pt_bits;

%% STBC Setup

rng('default');
rng('shuffle');

%QPSK symbols
st_pu_const = (exp(1i.*[-3*pi/4 3*pi/4 7*pi/4 -7*pi/4])); %Unity magnitude and thus unity power
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


%% Code

m = 8;

snr_dB = -5:5:30;
pl_st_pr_dB = 5;

	st_pr_snr = snr_dB(m) +pl_st_pr_dB;
	st_pr_snr = 10.^(st_pr_snr/10);
	st_pr_sigma = sqrt(1/st_pr_snr);
    
    
pt_symbols = randi([0 pt_M-1], 1,1)

pt_bit_sequence = de2bi(pt_symbols,pt_bits)

%At secondary transmitter

pt_st_y_decoded = pt_symbols;

pt_st_y_decoded_bits = de2bi(pt_st_y_decoded,pt_bits);

st_symbols = randi([1,4], 3, 1);

st_bit_sequence = de2bi(st_symbols-1, 2);

st_symbols = st_su_const(st_symbols);

c = st_pu_const(bi2de(reshape(pt_st_y_decoded_bits, 3, 2)) + ones(3,1))

s = st_symbols;

st_stbc_code = [c(1) c(2) c(3); -c(2)' c(1)' s(1)'; -c(3)' s(2)' c(1)'; s(3)' -c(3)' c(2)'];



%At primary receiver
h = (randn(3,1)+1i*randn(3,1))./sqrt(2) %Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.

%the following complex equivalent channel matrix is for the ORTHOGONAL 3TX STBC
H = [h(1) h(2) h(3); h(2)' -h(1)' 0; h(3)' 0 -h(1)'; 0 h(3)' -h(2)'];

%the following complex equivalent channel is for the Embedded Diversity Code
%for 3 TX
H_eqv = [h(1) h(2) h(3) 0 0 0; h(2)' -h(1)' 0 h(3)' 0 0; h(3)' 0 -h(1)' 0 h(2)' 0; 0 h(3)' -h(2)' 0 0 h(1)'];

st_pr_sigma
N = st_pr_sigma.*((randn(4,1)+1i*randn(4,1))./sqrt(2));

% T= 4, M_r = 1
%Joint variance of complex Gaussian distribution is 1. Therefore, average value of magnitude of fading channel is 1.

%received Signal
Y = st_stbc_code * h + N; %Code multiplied with h, not H i.e. multiplied with fading coefficients only.

for k = 1:length(X)
    %discard the effect of diversity 2 and 1 layer from total received
    %signal to get Y_remaining, alias Y_rem
    Y_rem = [Y(1); Y(2)-h(3)*X(k,1)'; Y(3)-h(2)*X(k,2)'; Y(4)-h(1)*X(k,3)'];
    
    Y_prime = [Y_rem(1) Y_rem(2:4)'].';
    
    %Apply matched filtering because the remaining received signal is
    %due to the contribution from the Orthogonal Diversity 3 layer Only
    Y_match = H' * Y_prime;
    
    Sym = sign([real(Y_match); imag(Y_match)]);
    
    S_tilde = Sym(1:3) + 1i*Sym(4:6);
    
    Decoded_Symb{k} = [S_tilde.' X(k,1) X(k,2) X(k,3)];
    %Now apply ML decoding using the overall equivalent channel matrix
    %H_Eqv
    diff = [Y(1) Y(2:4)'].' - H_eqv * Decoded_Symb{k}.';
    
    metric(k) = norm(diff,'fro')^2;
    
end

[W, ind] = min(metric);

st_pr_decoded_stbc = Decoded_Symb{ind};

st_pr_decoded_symbols = st_pr_decoded_stbc(1:3);

st_pr_decoded_bits = reshape([(sign(imag(st_pr_decoded_symbols))+1)/2;(sign(real(st_pr_decoded_symbols))+1)/2].', 1,6)
				
errors_pu_relay = sum(pt_bit_sequence ~= st_pr_decoded_bits)