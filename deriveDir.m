function[am,sm,dm,mod]=deriveDir(a,s,d)
load('data.mat');
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*exp(-((t-d).^2)/(2*(1+s*s).^2));
end
am=k;
k=0;
for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+(-4)*s/((1+s*s).^3)*((t-d).^2)*f*(sig_noisy(1,t)-f);
end
sm=k;
k=0;
for t=1:100
        k=k+(-2)*a*a*(exp(-(t-d)^2/(2*(1+s^2)^2)))...
        *(sig_noisy(1,t)-a^2*exp(-(t-d)^2/(2*(1+s^2)^2)))*(t-d)/((1+s^2)^2);
end
dm=k;
mod=sqrt(am*am+sm*sm+dm*dm);