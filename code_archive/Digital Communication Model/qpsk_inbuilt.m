%Information source
    sample_size = 100;
    bit_sequence = randi([0 1], [1 sample_size]);

%Source encoder: Express information source in minimum possible bits

    source_encoded_sequence = bit_sequence;

%Channel encoder: Binary data to code words with added error protection

    word_length = 2;

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
    
    modulator = comm.PSKModulator('ModulationOrder', constellation_elements, 'PhaseOffset',pi/4);
    
    modulated_baseband = step(modulator, symbols);
    % Obsololete function: pskmod(symbols,constellation_elements, pi/4);

    
%Channel

    EbNo = 10;
    
    awgn_noise_object = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', 'EbNo',  EbNo);
    recieved_signal = step(awgn_noise_object, modulated_baseband);
    
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
        
        
        %Reciever filter
            
        
        %Equalizer
        
        %Constellation diagram
            constellation_diagram = scatterplot(recieved_signal, 1, 0, 'g.');
            hold on
            scatterplot(modulated_baseband, 1, 0, 'k*', constellation_diagram);


    %Detection
    
        %Detection
        
        demodulator = comm.PSKDemodulator('ModulationOrder', constellation_elements, 'PhaseOffset',pi/4);
        detected_symbols = step(demodulator, recieved_symbols);
        
        %Obsolete demodulator: pskdemod(recieved_symbols, constellation_elements, pi/4);

        detected_words = de2bi(detected_symbols, word_length);

        detected_sequence = reshape(detected_words', [], 1)';

        recieved_sequence = detected_sequence;


        %BER Calculations
        [numErrors, ber] = biterr(source_encoded_sequence, recieved_sequence);
        fprintf('\nBER = %5.2e, %d errors', ber, numErrors);

%Channel Decoder

%Source Decoder