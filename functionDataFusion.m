function xhat = functionDataFusion(S,K,channelGaindB,numRealizations,x,dataFusion)
% Performs the data fusion step among the received signals at each subarray.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  S               number of subarrays.
% @param  K               number of users.
% @param  channelGaindB   M x K matrix with the channel gain in Decibels for each user when considering pathloss.
% @param  numRealizations number of channel realizations (small-fading).
% @param  x               K x 1 vector with modulated signals sent by all K users.
% @param  dataFusion      string that contains the data fusion type. Types are: 'DEDF' and 'DLDF'.
% @return xhat            K x numRealizations matrix with estimated signals at each channel realization.

%Prepare to save xhat
xhat = zeros(K,numRealizations);

%Check data fusion method
if strcmp(dataFusion,'DEDF')

    xhat = mean(x,3);

elseif strcmp(dataFusion,'DLDF')

    %Converting dB to linear scale
    channelGain = db2pow(channelGaindB);
    channelGain(channelGain == 1) = 0;

    %Store the sum of the channel gains
    overallSNR = reshape(sum(sum(channelGain,1),3),[K 1]);

    %Prepare to save the weigths
    alpha = zeros(K,S);

    %Go through each subarray
    for s = 1:S

        %Store the sum for each user
        subarraySNR = sum(channelGain(:,:,s),1);

        %Go through each user
        for k = 1:K

            %Obtain the weight by performing the ratio
            alpha(k,s) = subarraySNR(k)./overallSNR(k);

        end

    end

    %Go through each subarray
    for s = 1:S

        %Go through each user
        for k = 1:K

            xhat(k,:) = xhat(k,:) + alpha(k,s)*x(k,:,s);

        end

    end

end
