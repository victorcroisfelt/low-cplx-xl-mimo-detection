function xHat = functionSignalEstimate(M,K,numRealizations,Y,V,numIterRange)
% Outputs signal estimates. It is possible to enter with different numbers of iterations.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  M               number of antennas.
% @param  K               number of users.
% @param  numRealizations number of channel realizations (small-fading).
% @param  p               total uplink transmit power per UE [mW].
% @param  Y               M x numRealizations matrix with received signal.
% @param  V               M x numRealizations x K receive combining matrix.
% @param  numIterRange    vector with number of iterations points.
% @return xHat            matrix with estimate of the modulated signals sent by all K users; the estimate is based on the receive
% combining method.
%

%Check number of inputs
if nargin <= 5 % RZF

  %Prepare to save the signal estimates
  xHat = zeros(K,numRealizations);

  %Go through all channel realizations
  for n = 1:numRealizations

      %Extract current receive combining matrix
      Vn = reshape(V(:,n,:),[M K]);

      %Compute signal estimate
      xHat(:,n) = Vn'*Y(:,n);

  end

else % rKA-based

  %Prepare to save the signal estimates
  xHat = zeros(K,numRealizations,length(numInterRange));

  %Go through all channel realizations
  for n = 1:numRealizations

    %Go through all iterates
    for int = 1:length(numInterRange)

      %Extract current receive combining matrix
      Vn = reshape(V(:,n,:,int),[M K]);

      %Compute signal estimate
      xHat(:,n,int) = Vn'*Y(:,n);

    end

  end

end
