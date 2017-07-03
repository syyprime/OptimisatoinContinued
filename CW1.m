function[bool]= CW1(dk,xk,b1,a0)
load('data.mat');
gamma=b1*(dk*dk');
a=xk(1,1);
s=xk(1,2);
d=xk(1,3);
fxk=0;
for t=1:100
       fxk=fxk+(sig_noisy(1,t)-a*a*exp(-((t-d).^2)/(2*((1+s*s).^2)))).^2;
end
fxk1=0;
a1=a+a0*dk(1,1);
s1=s+a0*dk(1,2);
d1=d+a0*dk(1,3);
for t=1:100
       fxk1=fxk1+(sig_noisy(1,t)-a1*a1*exp(-((t-d1).^2)/(2*((1+s1*s1).^2)))).^2;
end
if(fxk1<=(fxk-a0*gamma))
    bool=1;
else
    bool=0;
end
