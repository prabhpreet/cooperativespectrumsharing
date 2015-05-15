sample_size = 100;
bit_rate = 1000; %bits per second
bit_rate_dB = 10*log10(bit_rate);
word_length = 1;
symbol_rate = bit_rate/word_length;
samples_per_symbol = 100; %Integer value, for simplicity
sampling_rate = samples_per_symbol*symbol_rate; %Hz
constellation_elements = 2.^word_length;

modulator = comm.PSKModulator('ModulationOrder', constellation_elements, 'PhaseOffset', 0);
demodulator = comm.PSKDemodulator('ModulationOrder', constellation_elements, 'PhaseOffset',0);

points_EbNo_dB = 15;
max_EbNo_dB = 15;
min_EbNo_dB = 0;
step_EbNo_dB = (max_EbNo_dB-min_EbNo_dB)/points_EbNo_dB;

runs = 100;

ber_awgn_only = zeros(runs,points_EbNo_dB);
ber_rayleigh = zeros(runs,points_EbNo_dB);

for j = 1:runs
	for i = 1:points_EbNo_dB
		EbNo_dB = min_EbNo_dB + (i*step_EbNo_dB);
		noise_power_dB = -EbNo_dB ;
		noise_power = 10^(noise_power/10);
		
		
	%Information source
		bit_sequence = randi([0 1], [1 sample_size]);

	%Source encoder: Express information source in minimum possible bits

		source_encoded_sequence = bit_sequence;

	%Channel encoder: Binary data to code words with added error protection

	   
		source_encoded_sequence = [source_encoded_sequence, zeros(1, mod(sample_size, word_length))]; %Pad zeros at end of odd sequence.

		words = reshape(source_encoded_sequence, word_length, [])';

		%Mapping scheme
		symbols = bi2de(words); %Binary mapping


	%Modulation

		modulated_baseband = step(modulator, symbols);
		modulated_baseband = modulated_baseband./abs(modulated_baseband);
		% Obsolete function: pskmod(symbols,constellation_elements);

		
	%Channel
		
		
		complex_fading_coefficient = (randn(length(symbols),1) + 1i.*rand(length(symbols),1))./sqrt(2);
		
		awgn_noise = noise_power.*(randn(length(symbols),1) + 1i.*rand(length(symbols),1))./sqrt(2);
		
		recieved_signal =  modulated_baseband + awgn_noise;
		
		awgn_noise =  noise_power.*(randn(length(symbols),1) + 1i.*rand(length(symbols),1))./sqrt(2);
		
		recieved_rayleigh_signal =  (complex_fading_coefficient.* modulated_baseband) + awgn_noise;
		
		
		%Adding AWGN: Obsolete
	% 
	%     snr_dB = EbNo_dB + 10*log10(word_length) - 10*log10(samples_per_symbol);
	% 
	%     recieved_signal = awgn(modulated_baseband, snr_dB, 'measured');

	%Demodulator

		%Demodulation
			
			%Frequency downconversion
				%Bandpass filter, etc.
			
			recieved_symbols = recieved_signal;
			recieved_rayleigh_symbols = recieved_rayleigh_signal;
			
			%Reciever filter
				
			
			%Equalizer
		
		%Detection
		
			%Detection

			
			detected_symbols = step(demodulator, recieved_symbols);
			detected_rayleigh_symbols = step(demodulator, recieved_rayleigh_symbols);
			
			%Obsolete demodulator: pskdemod(recieved_symbols, constellation_elements);

			detected_words = de2bi(detected_symbols, word_length);

			detected_sequence = reshape(detected_words', [], 1)';

			recieved_sequence = detected_sequence;

			detected_rayleigh_words = de2bi(detected_rayleigh_symbols, word_length);

			detected_rayleigh_sequence = reshape(detected_rayleigh_words', [], 1)';

			recieved_rayleigh_sequence = detected_rayleigh_sequence;
			
			%BER Calculations
			[numErrors, ber_awgn_only(j,i)] = biterr(source_encoded_sequence, recieved_sequence);

			
			[numRayleighErrors, ber_rayleigh(j,i)] = biterr(source_encoded_sequence, recieved_rayleigh_sequence);
	end
end

ber_awgn_average = sum(ber_awgn_only)./runs;

figure
EbNo_dB = min_EbNo_dB:step_EbNo_dB:max_EbNo_dB-step_EbNo_dB;
semilogy(EbNo_dB, ber_awgn_average, 'go-');


%hold on
%plot(EbNo_dB, ber_rayleigh, 'bo-');
%Channel Decoder

%Source Decoder