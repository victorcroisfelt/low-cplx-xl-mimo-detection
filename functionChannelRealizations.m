function H = functionChannelRealizations(M,K,channelGaindB,R,numRealizations)
% Generates channel realizations. Assumes Rayleigh fading and covariance matrix R.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of BS antennas.
% @param  K               number of users.
% @param  channelGaindB   M x K matrix with the channel gain in Decibels for each user when considering pathloss.
% @param  R               M x M x K matrix with K x (M x M) covariance matrices of K users.
% @param  numRealizations number of channel realizations (small-fading).
% @return H               M x numRealizations x K matrix with channel responses.
%

%% Channel realizations

%Generate small-scale fading realizations
Hbar = sqrt(0.5)*(randn(M,numRealizations,K)+1i*randn(M,numRealizations,K));

%Converting dB to linear scale
channelGain = db2pow(channelGaindB);
channelGain(channelGain == 1) = 0;

%Prepare to store channel matrix
H = zeros(M,numRealizations,K);

%Go through all users
for k = 1:K

    if trace(R(:,:,k)) ~= 0

        %Obtain the sqrt of channel covariance matrix
        Rsqrt = sqrtm(R(:,:,k));

        %Apply channel covariance structure to Hbar
        Hbar(:,:,k) = Rsqrt*Hbar(:,:,k);

        %Apply channel gain to Hbar and obtain H
        H(:,:,k) = sqrt(channelGain(:,k)).*Hbar(:,:,k);

    end

end
