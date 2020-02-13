function sumSINR = functionComputeSINR(M,K,p,numRealizations,H,V)
% Computes the sum of the signal-to-interference-plus-noise-ratio (SINR) of each user.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of BS antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  H               M x numRealizations x K matrix with channel responses.
% @return V               M x numRealizations x K RZF receive combining matrix.
%

%% Preamble

%Hold the K x K identity matrix
eyeK = eye(K);

%Prepare to save average sum SINR
sumSINR = zeros(K,numRealizations);

%Go through all channel realizations
for n = 1:numRealizations

    %Extract channel realizations from all users to the BS
    Hn = reshape(H(:,n,:),[M K]);

    if nargin <= 5

        %Compute RZF combining
        V_RZF = p*Hn/(p*(Hn'*Hn)+eyeK);

    end

    %Go through all users
    for k = 1:K

        %Extract receive combining vector
        if nargin <= 5

            v = V_RZF(:,k);

        else

            v = V(:,n,k);

        end

        %Compute signal and interference + noise terms
        signal = p*abs(v'*Hn(:,k)).^2;
        combiningNorm = norm(v).^2;
        interf = p*sum(abs(v'*Hn).^2);

        %Compute sum SINR
        sumSINR(k,n) = signal / (interf - signal + combiningNorm);

    end

end

%Treating NaNs
sumSINR(isnan(sumSINR)) = 0;

%Averaging and summing
sumSINR = sum(mean(sumSINR,2));
