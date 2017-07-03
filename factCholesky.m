function[L,tau]=factCholesky(hess)
hessT=hess;
H=norm(hess,'fro');
hii=min([hess(1,1),hess(2,2),hess(3,3)]);
if(hii>0)
    tau=0;
else

    tau=H;
end
flag=0;
lamda=eig(hessT);
for i=1:size(lamda)
    if(lamda(i,1)<0)
        flag=1;
    end
end

while(flag==1)
    flag=0;
    hessT=hess+eye(3).*tau;
    lamda=eig(hessT);
    for i=1:size(lamda)
        if(lamda(i,1)<0)
            flag=1;
        end
    end
    if(flag==0)
        break;
    else
        tau=max([2*tau,H/2]);
    end
end
L=chol(hessT);