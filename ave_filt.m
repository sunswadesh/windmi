% Backward average of n pts
function Y = ave_filt(X,n)

P = X(1,1)*ones(1,n-1);
FM = [P X];
filt_B = ones(1,n);
filt_A = [n];
Out = filter(filt_B,filt_A,FM);
Y = Out(n:size(FM,2));
