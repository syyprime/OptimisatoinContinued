function f=myfun(theda)

load('data.mat');
a=theda(1);
d=theda(2);
s=theda(3);
f=0;

for i=0:99
    f = f + (sig_noisy(i+1)-a.^2.*exp(-(i-d).^2./(2.*(1+s.^2).^2))).^2; 
end

end