function CRD = functionCRD(actualValue,upperBound)
% Computes computational relaxation factor (CRD).
% CRD = (upperBound - actualValue)/upperBound if actualValue is less than the upperBound; otherwise CRD = 0.
% <p>
% @author Victor Croisfelt <victorcroisfelt@gmail.com>
% </p>
% @param  actualValue     vector with current values.
% @param  upperBound      upper bound.
% @return CRD             vector of size of actualValue with computational relaxation factors.
%

%Get size of actual value
s = size(actualValue);

%Compute CRD
CRD = (repmat(upperBound,s) - actualValue)./repmat(upperBound,s);
CRD(CRD < 0) = 0;
