function [r, errbr] = calc_err(x,y)

cor=corrcoef(x,y);
r =cor(1,2);
fr= 0.5*log((1+r)/(1-r));
ufr= fr + 1.96*1/sqrt(length(y)-3);
lfr= fr - 1.96*1/sqrt(length(y)-3);
ur = (exp(2*ufr)-1)/(exp(2*ufr)+1);
lr = (exp(2*lfr)-1)/(exp(2*lfr)+1);
errbr=[ur-r r-lr];