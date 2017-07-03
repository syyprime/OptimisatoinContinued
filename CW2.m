function[bool]= CW2(dk,xk,b2,a0)
a=xk(1,1);
s=xk(1,2);
d=xk(1,3);
a1=a+a0*dk(1,1);
s1=s+a0*dk(1,2);
d1=d+a0*dk(1,3);
[fd1(1,1),fd1(1,2),fd1(1,3),mod1]=deriveDir(a1,s1,d1);
up=fd1*dk';
[fd(1,1),fd(1,2),fd(1,3),mod]=deriveDir(a,s,d);
down=fd*dk';
if(up/down<=b2)
    bool=1;
else
    bool=0;
end
