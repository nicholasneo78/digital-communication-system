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
AverageOOKError = zeros(1,11);

S=1;
error = 0;

for runs = 1:20
    bitErrorRateOutput = zeros(1,11);
    arrayIndex=1;
    for i = SNR
        SNRvalues(arrayIndex) = i;
        
        N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
        outputSignal = awgn(sampled_ook,i,N); %i=SNR, N=noise variance
        
        demodulated_signal = outputSignal .* (2 * carrier);
        
        filtered_signal = filtfilt(b, a, demodulated_signal);
        
        decoded_signal = zeros(1,1024);
        
        for count=1:noOfBits
            produced_signal = filtered_signal(1 /2 * Fs/dataRate + (count - 1) * Fs/dataRate);
            if (produced_signal > 0.5)
                decoded_signal(count) = 1;
            else
                decoded_signal(count) = 0;        
            end
            
            %think about how to put the error count into here from p1
            %instead of calling the function
        end
        
        SNRvalues(arrayIndex)=i; %Store current SNR value into array
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(arrayIndex)= bitErrorRate;
        AverageOOKError(arrayIndex) = AverageOOKError(arrayIndex) + bitErrorRateOutput(arrayIndex);
        arrayIndex = arrayIndex + 1;
    end
end

AverageOOKError = AverageOOKError ./ 20;

semilogy(SNRvalues, AverageOOKError);
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

for runs = 1:20
    bitErrorRateOutput = zeros(1,11);
    arrayIndex=1;
    for i = SNR
        SNRvalues(arrayIndex) = i;
        
        N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
        outputSignal = awgn(sampled_bpsk,i,N); %i=SNR, N=noise variance

        demodulated_signal = outputSignal .* (2 * carrier);
        
        filtered_signal = filtfilt(b, a, demodulated_signal);
        
        decoded_signal = zeros(1,1024);
        
        for count=1:noOfBits
            produced_signal = filtered_signal(1 /2 * Fs/dataRate + (count - 1) * Fs/dataRate);
            if produced_signal > 0
                decoded_signal(count) = 1;
            else
                decoded_signal(count) = 0;        
            end
        end
        
        SNRvalues(arrayIndex)=i; %Store current SNR value into array
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(arrayIndex)=bitErrorRate;
        AverageBPSKError(arrayIndex) = AverageBPSKError(arrayIndex) + bitErrorRateOutput(arrayIndex);
        arrayIndex = arrayIndex +1;
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