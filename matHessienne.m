function[hess]=matHessienne(a1,s1,d1,precision)
load('data.mat');
hess=zeros(3,3);
k=0;
a=a1;
s=s1;
d=d1;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*exp(-((t-d).^2)/(2*(1+s*s).^2));
end
am1=k;
%?f(¦Èk)/(?a*?a)
a=a1+precision;
s=s1;
d=d1;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*exp(-((t-d).^2)/(2*(1+s*s).^2));
end
ama=k;
hess(1,1)=(ama-am1)/precision;
%?f(¦Èk)/(?a*?s)
a=a1;
d=d1;
s=s1+precision;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*exp(-((t-d).^2)/(2*(1+s*s).^2));
end
ams=k;
hess(1,2)=(ams-am1)/precision;
%?f(¦Èk)/(?a*?d)
a=a1;
d=d1+precision;
s=s1;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*exp(-((t-d).^2)/(2*(1+s*s).^2));
end
amd=k;
hess(1,3)=(amd-am1)/precision;
%?f(¦Èk)/(?s*?a)
hess(2,1)=hess(1,2);
%?f(¦Èk)/(?s*?s)
a=a1;
d=d1;
s=s1;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+(-4)*s/((1+s*s).^3)*((t-d).^2)*f*(sig_noisy(1,t)-f);
end
sm=k;
a=a1;
d=d1;
s=s1+precision;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+(-4)*s/((1+s*s).^3)*((t-d).^2)*f*(sig_noisy(1,t)-f);
end
sms=k;
hess(2,2)=(sms-sm)/precision;
%?f(¦Èk)/(?s*?d)
a=a1;
d=d1+precision;
s=s1;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+(-4)*s/((1+s*s).^3)*((t-d).^2)*f*(sig_noisy(1,t)-f);
end
smd=k;
hess(2,3)=(smd-sm)/precision;
%?f(¦Èk)/(?d*?a)
hess(3,1)=hess(1,3);
%?f(¦Èk)/(?d*?s)
hess(3,2)=hess(2,3);
%?f(¦Èk)/(?d*?d)
a=a1;
s=s1;
d=d1;
k=0;
for t=1:100
        k=k+(-2)*a*a*(exp(-(t-d)^2/(2*(1+s^2)^2)))...
        *(sig_noisy(1,t)-a^2*exp(-(t-d)^2/(2*(1+s^2)^2)))*(t-d)/((1+s^2)^2);
end
dm=k;
a=a1;
s=s1;
d=d1+precision;
k=0;
for t=1:100
        k=k+(-2)*a*a*(exp(-(t-d)^2/(2*(1+s^2)^2)))...
        *(sig_noisy(1,t)-a^2*exp(-(t-d)^2/(2*(1+s^2)^2)))*(t-d)/((1+s^2)^2);
end
dmd=k;
hess(3,3)=(dmd-dm)/precision;

