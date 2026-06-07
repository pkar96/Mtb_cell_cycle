function [cv, errbr] = calc_cv(x)

cv= std(x)/mean(x);

errbr=cv^2*sqrt((1+sqrt(2)*(mean(x)/std(x))^2)/length(x));