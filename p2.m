%Define carrier frequency in Hz
Fc = 10000; 
%Given 16 times over-sampled Fs = Fc * 16
Fs = Fc * 16;
%Data rate 1kpbs
dataRate = 1000;
%Signal length
noOfBits = 1024;
%sampling time
sampStart = 1/(2 * Fs);
sampInterval = 1/Fs;
timeTaken = noOfBits/dataRate;
time = sampStart: sampInterval: timeTaken;
%carrier
carrier = cos(2 * pi * Fc * time);

%input
input = randi([0, 1], [1, noOfBits]);
%Simulation for different SNR values
SNR = 0:5:50; %Loop from 0 to 50 (in multiples of 5)

%Define amp
amplitude = 1;

%low pass butter filter, 6th order filter with cut-off freq 0.2
[b, a] = butter(6, 0.2);



%sampling rate is larger than data rate
%need to extend 1s and 0s by the ratio of the sampling rate and the data
%rate
ratio_fs_dataRate = Fs/dataRate;
extension = ones(1, ratio_fs_dataRate);
sampled_input = kron(input, extension);

%extended sampled input multipled with carrier signal
sampled_ook = sampled_input .* carrier;

SNRvalues = zeros(1,11); %Initialise an array of 1-by-11 zeros to store SNR values
AveragebitErrorRateOutput = zeros(1,11);

for runs = 1:20
    bitErrorRateOutput = zeros(1,11);
    counter=1;
    for SNR = 0:5:50
        SNRvalues(counter) = SNR;
        noise_variance = 1 / 10^(SNR/10);
        noise_std = sqrt(noise_variance);
        noise = noise_std .* randn(1, 1024 * Fs/dataRate);
        
        received_signal = sampled_ook + noise;
        demodulated_signal = received_signal .* (2 * carrier);
        filtered_signal = filtfilt(b, a, demodulated_signal);
        decoded_signal = zeros(1,1024);
        for i = 1:1:1024
            interested_signal = filtered_signal(1 /2 * Fs/dataRate + (i - 1) * Fs/dataRate);
            if interested_signal > 0.5
                decoded_signal(i) = 1;
            else
                decoded_signal(i) = 0;        
            end
        end
        
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(counter)= bitErrorRate;
        AveragebitErrorRateOutput(counter) = AveragebitErrorRateOutput(counter) + bitErrorRateOutput(counter);
        counter = counter +1;
    end
end

AveragebitErrorRateOutput = AveragebitErrorRateOutput ./ 20;

semilogy(SNRvalues, AveragebitErrorRateOutput);
ylim([10^(-5) 10^1]);
xlim([0 50]);
hold on

hold on
S=1;
SNRvalues = zeros(1,11);
bitErrorRateOutput = zeros(1,11);
counter=1;
AverageBPSKError = zeros(1,11);

%test
% inputData = randi([0, 1], [1, 2048]);
% ratio_fs_dataRate = Fs/dataRate;
% extension = ones(1, ratio_fs_dataRate);
% sampled_input = kron(inputData, extension);

%extended sampled input multipled with carrier signal
sampled_input_bpsk = 2 * sampled_input - 1;
sampled_bpsk = sampled_input_bpsk .* carrier;
arrayIndex = 1; %To count the array index for storage

for i=SNR
    noOfErrors = 0; %Count number of errors for each SNR
    
    N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
    
    %AWGN function applies White Gaussian Noise to inputData
    outputSignal = awgn(sampled_bpsk,i,N); %i=SNR, N=noise variance
    
    %Convert outputSignal to 0/1
    %(outputSignal >= 0) -> 1
    %(outputSignal < 0) -> 0
    threshold = 0; 
    outputSignal(outputSignal>=threshold) = 1;
    outputSignal(outputSignal<threshold) = 0;
end
for runs = 1:20
    bitErrorRateOutput = zeros(1,11);
    counter=1;
    for SNR = 0:5:50
        SNRvalues(counter) = SNR;
        noise_variance = 1 / 10^(SNR/10);
        noise_std = sqrt(noise_variance);
        noise = noise_std .* randn(1, 1024 * Fs/dataRate);
        
        received_signal = sampled_bpsk + noise;
        demodulated_signal = received_signal .* (2 * carrier);
        filtered_signal = filtfilt(b, a, demodulated_signal);
        decoded_signal = zeros(1,1024);
        for i = 1:1:1024
            interested_signal = filtered_signal(1 /2 * Fs/dataRate + (i - 1) * Fs/dataRate);
            if interested_signal > 0
                decoded_signal(i) = 1;
            else
                decoded_signal(i) = 0;        
            end
        end
        
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(counter)=bitErrorRate;
        AverageBPSKError(counter) = AverageBPSKError(counter) + bitErrorRateOutput(counter);
        counter = counter +1;
    end
end
AverageBPSKError = AverageBPSKError ./ 20;
semilogy(SNRvalues, AverageBPSKError);
ylabel('Log 10 Bit Error Rate') ;
hold on
title('Plot of Bit Error vs SNR for OOK and BPSK');
legend({'y = AverageOOK','y= AverageBPSK'},'Location','southeast')
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;

function bitErrorRate = calculate_error_rate(input, tempInput)
    %Generate noise having normal distribution with zero mean
    error = 0;
    numOfBits= 1024;
    for i = 1:1:numOfBits
        if input(i) ~= tempInput(i)
            error = error + 1;
        end
    end
    bitErrorRate = error/numOfBits;
end