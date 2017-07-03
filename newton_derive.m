function d_fin = newton_derive(f,f1,f2,d_debut)
syms d;int i=1;int cout_min=0;
c=zeros(1,100);
f_cout(d)=f;
f_derive(d)=f1;
f_2derive(d)=f2;
while(d_debut<=100)
    d_fin=d_debut-double(f_derive(d_debut))/double(f_2derive(d_debut));
    d_debut=d_debut+1
    if(double(f_derive(d_fin))>=1e-9 && f_2derive(d_fin)>0)
        c(i)=d_fin;
        i=i+1;
    end
end
    
for x=1:100
    di=c(x);
    cout_min=double(f_cout(di));
    if(cout_min>double(f_cout(di)))
        cout_min=double(f_cout(di));
        d_fin=c(x);
    end
end
return