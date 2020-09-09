%[1-1]
% generate random binary digits (0 or 1) 
noOfBits = 1024;
inputData = randi([0, 1], [1, noOfBits]);
%inputData % test data generated

%[1-2]
% convert binary digits to +/- 1 (1->+1, 0->-1) => Data for transmission
% replace 0 to -1
inputData(inputData==0) = -1;
%inputData % test data generated

zeroMean = 0;
SNR = 0:5:50; %Loop from 0 to 50 (in multiples of 5)
S=1;
bitErrorRate = zeros(1,11); %Initialise an array of 1-by-11 zeros
SNRvalues = zeros(1,11);
arrayIndex = 1;

for i = SNR
    noOfErrors = 0;
    N=S./(10.^(i./10));
    %noise = wgn(1,1024,sqrt(N));
    %outputSignal = awgn(inputData,i,N,(randn(1024,1)*sqrt(1)));
    %calculate = sum(noise);
    %calculate
    %outputSignal = noise+inputData;
    threshold = 0;
    outputSignal(outputSignal>=0) = 1;
    outputSignal(outputSignal<0) = 0;
    plot(outputSignal);
    for count=1:noOfBits
        if (inputData(count) > threshold && outputSignal(count) == 0) || (inputData(count) <= threshold && outputSignal(count) == 1)
            noOfErrors = noOfErrors + 1;
        end
    end
    SNRvalues(arrayIndex)=i;
    bitErrorRate(arrayIndex)= noOfErrors/noOfBits; %Store bitErrorRate into array
    arrayIndex=arrayIndex+1;
end

%Plot
semilogy(SNRvalues, bitErrorRate);
axis([0 50 -1 1])
title("Plot of Bit Error Rate vs Signal to Noise Ratio");
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;

%{

%[1-3, 1-4]
% generate equal number of noise samples
S = 1; %Assume signal (input data) has unit power
meanNoise = 0; % zero mean
SNR = 0:5:50; %Loop from 0 to 50 (in multiples of 5)
bitErrorRateOutput = zeros(1,11); %Initialise an array of 1-by-11 zeros
SNRAxis = zeros(1,11); %Initialise an array of 1-by-11 zeros
counter = 1; %Counter for array index

for i=SNR 
    varianceNoise = S./(10.^(i./10)); % unit variance
    % Follow notes
    receivedSignal = calculateNoise(inputData, varianceNoise);
    threshold = 0;
    receivedSignal(receivedSignal>=0) = 1;
    receivedSignal(receivedSignal<0) = 0;
    bitError = 0;
    for count=1 : noOfBits
        if (inputData(count) > threshold && receivedSignal(count) == 0) || (inputData(count) <= threshold && receivedSignal(count) == 1)
            bitError = bitError + 1;
        end
    end    
    bitErrorRateOutput (counter)= bitError/noOfBits; %Store bitErrorRate into array
    SNRAxis(counter) = i; %Store SNR value
    counter=counter+1;
end

%Plot
semilogy(SNRAxis, bitErrorRateOutput);
axis([0 50 -1 1])
title("Plot of Bit Error Rate vs Signal to Noise Ratio");
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;

function receivedSignal = calculateNoise(inputData, varianceNoise)
    % get noise (according to normal distribution)
    noise = (sqrt(varianceNoise).*randn(1,1024)); %Generate white gaussian noise
    %noise % look at the noise values generated with mean=0 and variance=1
    receivedSignal = noise+inputData; %Adding generated noise with generated random data
end

%}