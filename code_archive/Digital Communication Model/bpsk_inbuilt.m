%Information source
    sample_size = 100;
    bit_sequence = randi([0 1], [1 sample_size]);

%Source encoder: Express information source in minimum possible bits

    source_encoded_sequence = bit_sequence;

%Channel encoder: Binary data to code words with added error protection

    word_length = 1;

    constellation_elements = 2.^word_length;

    source_encoded_sequence = [source_encoded_sequence, zeros(1, mod(sample_size, word_length))]; %Pad zeros at end of odd sequence.

    words = reshape(source_encoded_sequence, word_length, [])';

    %Mapping scheme
    symbols = bi2de(words); %Binary mapping


%Modulation

    bit_rate = 1000; %bits per second
    symbol_rate = bit_rate/word_length;
    samples_per_symbol = 100; %Integer value, for simplicity
    sampling_rate = samples_per_symbol*symbol_rate; %Hz

    %Pulse signal
    
    modulator = comm.PSKModulator('ModulationOrder', constellation_elements, 'PhaseOffset', 0);
    
    modulated_baseband = step(modulator, symbols);
    % Obsololete function: pskmod(symbols,constellation_elements);

    
%Channel
    
     EbNo = 10;
    
    rayleigh_channel_object = rayleighchan(1/symbol_rate, 0); %No doppler shift, flat rayleigh channel
    rayleigh_modulated_baseband = filter(rayleigh_channel_object, modulated_baseband);
    
    
     
    awgn_noise_object = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', 'EbNo',  EbNo);
    
    recieved_signal = step(awgn_noise_object, modulated_baseband);
    recieved_rayleigh_signal = step(awgn_noise_object, rayleigh_modulated_baseband);
    
    
    %Adding AWGN: Obsolete
% 
%     snr = EbNo + 10*log10(word_length) - 10*log10(samples_per_symbol);
% 
%     recieved_signal = awgn(modulated_baseband, snr, 'measured');

%Demodulator

    %Demodulation
        
        %Frequency downconversion
            %Bandpass filter, etc.
        
        recieved_symbols = recieved_signal;
        recieved_rayleigh_symbols = recieved_rayleigh_signal;
        
        %Reciever filter
            
        
        %Equalizer
        
        %Constellation diagram
            constellation_diagram = scatterplot(recieved_symbols, 1, 0, 'g.');
            hold on
            scatterplot(modulated_baseband, 1, 0, 'k*', constellation_diagram);

            constellation_rayleigh_diagram = scatterplot(recieved_rayleigh_symbols, 1, 0, 'g.');
            hold on
            scatterplot(modulated_baseband, 1, 0, 'k*', constellation_rayleigh_diagram);

    %Detection
    
        %Detection

        demodulator = comm.PSKDemodulator('ModulationOrder', constellation_elements, 'PhaseOffset',0);
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
        [numErrors, ber] = biterr(source_encoded_sequence, recieved_sequence);
        fprintf('\nBER = %5.2e, %d errors', ber, numErrors);
        
        
        [numRayleighErrors, berRayleigh] = biterr(source_encoded_sequence, recieved_rayleigh_sequence);
        fprintf('\nRayleigh BER = %5.2e, %d errors', berRayleigh, numRayleighErrors);
%Channel Decoder

%Source Decoder