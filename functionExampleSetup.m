function [channelGaindB,R,meanUsers] = functionExampleSetup(M,S,K,diagNorm)
% Generates mobile communication setup. Assumes a sigle cell with square shape being supplied by a base station (BS) equipped with M antennas.
% The signal processing is divided into S subarrays, each containing Ms = M/S antennas. The cell serves K single-antenna users.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of BS antennas.
% @param  S               number of subarrays.
% @param  K               number of users.
% @param  diagNorm        string that contains different normalizations of subarrays' covariance matrices. Choices are: 'Norm1' and 'Norm2'.
% @return channelGaindB   M x K matrix with the channel gain in Decibels for each user when considering pathloss.
% @return R               M x M x K matrix with K x (M x M) covariance matrices of K users.
% @return meanUsers       scalar with the average number of users per subarray.
%

%% System parameters

%Number of antennas per fixed-size subarray
Ms = M/S;

%% Pathloss model parameters

%Atenuation coefficient
Omega = 4;

%Pathloss exponent
nu = 3;

%% Geometric and operational parameters

%Set the length in meters of the total square area
squareLength = 100;

%Operation frequency [Hz]
freq = 2.6e9;

%Operation Wavelength [m]
lambda = physconst('LightSpeed')/freq;

%Define the antenna spacing as half the wavelength
antennaSpacing = 2*lambda;

%Array size
arrayLength = M*antennaSpacing;

%Antenna positions
antennaPositions = [flip(-antennaSpacing:-antennaSpacing:-arrayLength/2) antennaSpacing:antennaSpacing:arrayLength/2];

%% Visibility regions parameters

% Mean of VR length
muVR = 0.1*arrayLength; %m
sigmaVR = 0.1; %m

%% Randomly putting out users in the cell

%Prepare to put out users in the cell
UEpositions = zeros(K,1);

%Auxiliar variable
perBS = 0;

%Put out K users in the cell, uniformly at random. The procedure is
%iterative since users that do not satisfy the minimum distance are
%replaced with new users
while perBS<K

    %Put out new users
    UEremaining = K - perBS;
    posX = rand(UEremaining,1)*squareLength - squareLength/2;
    posY = rand(UEremaining,1)*squareLength - squareLength/2;
    posXY = posX + 1i*posY;

    %Keep those that satisfy the minimum distance
    posXY = posXY(abs(posXY) >= 30);

    %Store new users
    UEpositions(perBS+1:perBS+length(posXY)) = posXY;
    perBS = perBS + length(posXY);

end

%% Generate mobile communication setup

%Generating VR centers
cVR = arrayLength*rand(K,1);

%Generating VR lengths
lVR = lognrnd(muVR,sigmaVR,[K 1]);

%Average number of users per subarray
meanUsers = 0;

%Prepare to save channel gains of each user
channelGaindB = zeros(Ms,K,S);

%Diagonal matrix denoting spatial non-stationarities
D = zeros(Ms,Ms,K,S);

%Prepare to store channel covariance matrices
R = zeros(Ms,Ms,K,S);

%Go through all users
for k = 1:K

    %Calculate distances between user k and BS antenna m
    distances = abs(UEpositions(k) - antennaPositions);

    %Obtain channel gain for the kth user
    channelGaindB(:,k,:) = reshape(-10*log10(Omega*distances.^nu),[Ms S]);

    %Search for antennas in the range of the VR
    antennaIndexes = ((cVR(k)-lVR(k) <= (antennaPositions+arrayLength/2)) + ((antennaPositions+arrayLength/2) <= cVR(k)+lVR(k))) > 1;

    if strcmp(diagNorm,'Norm1')

        %Normalization 1: trace of Omega must be equal to M (see paper)
        antennaIndexes = (M/sum(antennaIndexes))*antennaIndexes;

    end

    %Reshaping
    antennaIndexes = reshape(antennaIndexes,[Ms S]);

    %Go through each subarray
    for s = 1:S

        %Cleaning unusual channel gains
        channelGaindB((antennaIndexes(:,s) == 0),k,s) = 0;

        %Update average number of users per subarray
        meanUsers = meanUsers + double(sum(antennaIndexes(:,s)) > 0)/S;

        %Store diagonal matrix
        D(:,:,k,s) = diag(antennaIndexes(:,s));

        %Set the covariance matrix for each subarray
        R(:,:,k,s) = D(:,:,k,s).^(1/2)*eye(Ms)*D(:,:,k,s).^(1/2);

    end

end
