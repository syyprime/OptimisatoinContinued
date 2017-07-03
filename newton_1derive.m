function d_fin = newton_1derive(f1,f2,d0)
syms d;
f_derive(d)=f1;
f_2derive(d)=f2;
d_debut=d0;
i=1;
for i=1:100
    d_fin=d_debut-double(f_derive(d_debut))/double(f_2derive(d_debut));
    if(abs(d_fin-d_debut)>1)
        d_debut=d_fin;
    else
        break
    end
end
return