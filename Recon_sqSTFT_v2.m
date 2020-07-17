function [recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, Hz, c, Band) 

% Band = absolute frequency band usually is chosen to be 0.1Hz
% tfrsqtic = in the scale of 1Hz
% c is the extracted curve.

coeff=0.1152; % computed in the numerical example 'genFig1Fig2_test_coeff_SST.m'

alpha = tfrsqtic(2)-tfrsqtic(1) ;
RR = round(Band/(Hz*alpha));

recon = [] ;

for kk = 1: length(c)
	idx = max(1,c(kk)-RR): min(length(tfrsqtic),c(kk)+RR) ;
	recon(kk) = 2*sum(tfrsq(idx,kk),1)*(Hz/2/0.5)*(alpha)/Hz ;%
end

recon=recon*(coeff)*Hz/2;
