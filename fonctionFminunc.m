function[fk]=fonctionFminunc(r)
load('data.mat');
a=r(1);
s=r(2);
d=r(3);
fk=0;
for t=1:100
       fk=fk+(sig_noisy(1,t)-a*a*exp(-((t-d).^2)/(2*((1+s*s).^2)))).^2;
end
