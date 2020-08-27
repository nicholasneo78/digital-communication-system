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

%[1-3]
% generate equal number of noise samples
meanNoise = 0; % zero mean
varianceNoise = 1; % unit variance
% get noise (according to normal distribution)
noise = sqrt(varianceNoise)*randn(1,noOfBits) + meanNoise;
%noise % look at the noise values generated with mean=0 and variance=1
