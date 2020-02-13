function Y = functionReceivedSignal(M,K,p,numRealizations,x,H)
% Outputs the received signal at subarray s.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of antennas.
% @param  K               number of users.
% @param  p               total uplink transmit power per UE [mW].
% @param  numRealizations number of channel realizations (small-fading).
% @param  x               K x 1 vector with modulated signals sent by all K users.
% @param  H               M x numRealizations x K matrix with channel responses.
% @return Y               M x numRealizations matrix with received signal.
%

%Prepare to save received signal
Y = zeros(M,numRealizations);

%Go through all channel realizations
for n = 1:numRealizations

    %Extract current channel matrix
    Hn = reshape(H(:,n,:),[M K]);

    %Compute received signal
    Y(:,n) = sqrt(p)*Hn*x;

end

%Generate normalized random Gaussian noise
noise = sqrt(0.5)*(randn(M,numRealizations)+1i*randn(M,numRealizations));

%Additive Gaussian noise
Y = Y + noise
