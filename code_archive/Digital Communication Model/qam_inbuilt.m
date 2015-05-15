%Information source
    sample_size = 1000000;
    bit_sequence = randi([0 1], [1 sample_size]);

%Source encoder: Express information source in minimum possible bits

    source_encoded_sequence = bit_sequence;

%Channel encoder: Binary data to code words with added error protection

    word_length = 4;

    constellation_elements = 2.^word_length;

    source_encoded_sequence = [source_encoded_sequence, zeros(1, mod(sample_size, word_length))]; %Pad zeros at end of odd sequence.

    words = reshape(source_encoded_sequence, word_length, [])';

    symbols = bi2de(words);

%Digital Modulator

modulated_baseband = qammod(symbols, constellation_elements, 0); %Binary coding


%Channel

%Adding AWGN
EbNo = 20;
number_of_samples_per_symbol = 1;

snr = EbNo + 10*log10(word_length) - 10*log10(number_of_samples_per_symbol)

recieved_signal = awgn(modulated_baseband, snr, 'measured');

%Constellation diagram
constellation_diagram = scatterplot(recieved_signal, 1, 0, 'g.');
hold on
scatterplot(modulated_baseband, 1, 0, 'k*', constellation_diagram);

%Reciever

%Demodulation

recieved_symbols = recieved_signal;

%Detection

detected_symbols = qamdemod(recieved_symbols, constellation_elements);

detected_words = de2bi(detected_symbols, word_length);

detected_sequence = reshape(detected_words', [], 1)';

recieved_sequence = detected_sequence;


%BER Calculations
[numErrors, ber] = biterr(source_encoded_sequence, recieved_sequence);
fprintf('\nBER = %5.2e, %d errors', ber, numErrors);


