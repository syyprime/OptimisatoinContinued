close all;clear all;
load('data.mat');
%% Q2 la fonction de cout C(d) pour d variant entre 0 et 100.
syms d;
f(d)=d;
% on ecrit une fonction f(d) qui contient une valeur initiale.
% cette caleur va etre supprim¨¦e a la fin.
a=2;s=2;
cout=zeros(1,100);

for i=0:99
    f(d) = f(d) + (sig_noisy(i+1)-a^2*exp(-(i-d)^2/(2*(1+s^2)^2)))^2;
    %y = 4*gaussmf(i,[5 d]);
    %cout(1,d) = cout(1,d) + (sig_noisy(1,i)-y)^2;
end
f(d)=f(d)-d;
% maintenant f(d) est la somme pour tout les points i dans fonction cout.

for d1=0:99
    cout(d1+1) = f(d1);
end
% on met differentes valeurs de d dans la fontion.
% la valeur de cout(d) indique la distance carr¨¦ quand la variable=d.
figure(1);stem(cout);xlabel('d');ylabel('amplitude');
title('Fonction Cout');

%% Q4 la d¨¦riv¨¦e de la fonction de co?t en comparant avec diff 
cout_d1_manuel=zeros(1,100);
syms d;
f1m(d)=d;
for i=0:99
    f1m(d) = f1m(d) + (-2)*a*a*(exp(-(i-d)^2/(2*(1+s^2)^2)))*(sig_noisy(i+1)-a^2*exp(-(i-d)^2/(2*(1+s^2)^2)))*(i-d)/((1+s^2)^2);
end
f1m(d)=f1m(d)-d;
for d1=0:99
    cout_d1_manuel(d1+1) = f1m(d1);
end
figure(2);stem(cout_d1_manuel);xlabel('d');ylabel('amplitude');
title('D¨¦riv¨¦ manuel de Fonction Cout');

f1(d)=diff(f(d),d);
%f1 est la d¨¦riv¨¦ de la fonction cout 
cout_d1=zeros(1,100);
%on cr¨¦e une matrice pour d¨¦poser les valeurs de diff¨¦rents d.
for d1= 0:99
    cout_d1(d1+1)=f1(d1);
end
figure(3);stem(cout_d1);xlabel('d');ylabel('amplitude');
title('D¨¦riv¨¦ de Fonction Cout');
% cout_d1 contient tous les valeurs de different d.
%% Q6 la d¨¦riv¨¦e seconde de la fonction de co?t en comparant avec diff 
cout_d2_manuel=zeros(1,100);
syms d;
f2m(d)=d;
for i=0:99
    f2m(d) = f2m(d) + (-2)*(a^2)*exp(-((i - d)^2/(2 * (1 + s^2)^2)))*(-(sig_noisy(i+1)-(a^2)*exp(-((i - d)^2/(2 * (1 + s^2)^2))))/((1+s^2)^2)+(sig_noisy(i+1)-2*(a^2)*exp(-((i - d)^2/(2 * (1 + s^2)^2))))*((i-d)/((1+s^2)^2))^2);
end
f2m(d)=f2m(d)-d;
for d1=0:99
    cout_d2_manuel(d1+1) = f2m(d1);
end
figure(4);stem(cout_d2_manuel);xlabel('d');ylabel('amplitude');
title('D¨¦riv¨¦ secondaire manuel de Fonction Cout');

f2(d)=diff(f1(d),d);
cout_d2=zeros(1,100);

for d1= 0:99
    cout_d2(d1+1)=f2(d1);
end
figure(5);stem(cout_d2);xlabel('d');ylabel('amplitude');
title('D¨¦riv¨¦ secondaire de Fonction Cout');

%% Q7 l¡¯approximation quadratique par le d¨¦veloppement de Taylor
x=0;
f_tl(d) = f(d-x) + f1(d-x)*x + x*x/2*f2(d-x);
taylor=zeros(1,100);
n = randi(99);
x=-1;dn=n;
for d1= n:-1:1
    taylor(d1)=f_tl(dn);
    dn=dn-1;
end
x=1;dn=n;
for d1= (n+1):100
    taylor(d1)=f_tl(dn);
    dn=dn+1;
end
figure(6);stem(taylor);xlabel('d');ylabel('amplitude');title('Taylor');

%% Q8 les it¨¦r¨¦s de l¡¯algorithme de Newton pour diff¨¦rents points de d¨¦part
% pour les fonctions dont la valeur ne peut pas etre 0,on peut seulement
% calculer racine pour les minimiseurs.
% pour les fonctions qui n'ont pas de racine,cette methode ne peut pas nous 
% aider trouver tous les minimus/maximums globals d¨¦termin¨¦s. 
% on s'arr¨ºte ¨¤ un point dont le d¨¦riv¨¦ est presque 0.
d_fin=newton(f(d),f1(d),n);
d_ite=zeros(1,100);
for i= 0:99
    d_ite(i+1)=newton(f1(d),f2(d),i);
end
figure(7);stem(d_ite);xlabel('d');ylabel('amplitude');
title('les racines requises par differents points de depart');

%% Q9 la valeur qu'on obtient en s¡¯initialisant o¨´ le signal est maximum
[max_sig,position]=max(sig_noisy);
position=position-1;
% le maximum du signal est ¨¤ 25 colonne,mais le fonction depart ¨¤ 0,donc on
% le moins 1.
point=newton(f1(d),f2(d),position);
point=round(point);
valeur=cout(point);
i(1:100) = 0:99;
courbe_gauss = a^2*exp(-(i-point).^2/(2*(1+s^2)^2));
figure(8);plot(i,courbe_gauss);hold on;plot(i,sig_noisy);



%% Q10
clc;clear all;close all;load('data.mat');

% repr¨¦sentation de la fonction de cout ¨¤ plusieurs param¨¨tres
[a,d,s]=meshgrid(1:1:10,1:5:100,20:1:40);% on initialise d'abord les diff¨¦rents param¨¨tres
v=0;
for i=0:99
    v=v+(sig_noisy(i+1)-a.^2.*exp(-(i-d).^2./(2.*(1+s.^2).^2))).^2; 
    % ¨¦criture de la fonction de cout
end
aslice = [3,5,7]; 
dslice = [25,50,75]; 
sslice = [25,30,35];
figure(9);
slice(a,d,s,v,aslice,dslice,sslice);colormap cool;
% dessin de la fonction pour les valeurs indiqu¨¦es pr¨¦c¨¦demment


c=zeros(4,39,49);
% initialisation de la matric C qui va servir ¨¤ stocker les valeurs de la
% fonction de cout
for a=1:1:4 
    for s=1:2:49
        for d=1:2:39
            v=0;
            for i=0:99
                v=v+(sig_noisy(i+1)-a.^2.*exp(-(i-d).^2./(2.*(1+s.^2).^2))).^2; 
            end
            c(a,d,s)=v;
            % les valeurs de la focntion de cout sont stock¨¦es dans la matrice c 
            % ¨¤ la position correspondante.
        end
    end
end

% on fixe a ,d puis s et on trace la fonction de cout correspondante
figure(10);
va=0;
[d_a,s_a] = meshgrid(10:2:39,0:2:49);
for i=0:99   
    va=va+(sig_noisy(i+1)-5.^2.*exp(-(i-d_a).^2./(2.*(1+s_a.^2).^2))).^2;
end
subplot(2,2,1);surf(d_a,s_a,va);title('a=5');
vd=0;
[a_d,s_d] = meshgrid(0:1:4,0:2:49);
for i=0:99   
    vd=vd+(sig_noisy(i+1)-a_d.^2.*exp(-(i-20).^2./(2.*(1+s_d.^2).^2))).^2;
end
subplot(2,2,2);surf(a_d,s_d,vd);title('d=20');
vs=0;
[a_s,d_s] = meshgrid(0:1:4,10:2:39);
for i=0:99   
    vs=vs+(sig_noisy(i+1)-a_s.^2.*exp(-(i-d_s).^2./(2.*(1+2^2).^2))).^2;
end
subplot(2,2,3);surf(a_s,d_s,vs);title('s=2');


%% Q11
% On cherche le minimum des valeurs de la fonction de cout calcul¨¦es pour le
% triplet( a,d,s)
cout_min=min(min(min(c(1:4,1:2:39,1:2:49))));
% il suffit de trouver le mininmum de la fonction de cout dans la matrice C 
% qui a stock¨¦ les valeurs correspandentes.
[a_min,d_min,s_min]=find(c==cout_min);
% on cherche ici le triplet (a,d,s) optimal correspondant ¨¤ la valeur de la
% fonction de cout trouv¨¦e juste avant.

%% Q12
% 3^3=27 fois calculs pour (d,a,s)
% donc si on avait n variables inconnues nous auront 3^n calculs
% s'il existe n param¨¨tre, il existe aussi 3^n valeurs de la fonction
% ensemble ¨¤ comparer pour trouver le minimun parmi les valeurs.
% Donc si le nombre n de param¨¨tres est inconnu,il doit 3^n calculs.Quant
% au triplet ,il doit 3^3=27 calculs.
                
%% Q14
% le d¨¦riv¨¦e par rapport ¨¤ "a"
% par la fontion "diff"
s=2;
d=25;
i=1;
for a=0:0.01:3
    k=0;
    for t=1:100
        k=k+(sig_noisy(1,t)-a^2*exp(-((t-d).^2)/50)).^2;
    end
    f2(1,i)=k;
    i=i+1;
end
h=0.01;
x=0:h:3;
Y=diff(f2)/h;
% par le facon manuel
j=1;
for a=0:0.01:3
    k=0;
    for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+2*(sig_noisy(1,t)-f)*(-2*a)*f/(a^2);
    end
    f2m(1,j)=k;
    j=j+1;
end
figure(11);

plot(x(:,1:length(Y)),Y,'c*',x,f2m,'k.');
legend('la fonction diff','la fonction manuelle');
title('par rapport ¨¤''a''');

% le d¨¦riv¨¦e par rapport ¨¤ "d"
% par la fontion "diff"
s=2;a=2;i=1;
for d=0:0.1:100
    k=0;
    for t=1:100
        k=k+(sig_noisy(1,t)-4*exp(-((t-d).^2)/50)).^2;
    end
    f1(1,i)=k;
    i=i+1;
end
h=0.1;
x=0:h:100;
Y=diff(f1)/h;
% par le facon manuel
j=1;
for d=0:0.1:100
    k=0;
    for t=1:100
        k=k+(-2)*a*a*(exp(-(t-d)^2/(2*(1+s^2)^2)))...
        *(sig_noisy(1,t)-a^2*exp(-(t-d)^2/(2*(1+s^2)^2)))*(t-d)/((1+s^2)^2);
    end
    f1m(1,j)=k;
    j=j+1;
end
figure(12);
plot(x(:,1:length(Y)),Y,'c*',x,f1m,'k.');
legend('la fonction diff','la fonction manuelle');
title('par rapport ¨¤ ''d''');


% le d¨¦riv¨¦e par rapport ¨¤ "s"
% par la fontion "diff"
a=2;
d=25;
i=1;
for s=0:0.01:5
    k=0;
    for t=1:100
        k=k+(sig_noisy(1,t)-a^2*exp(-((t-d).^2)/(2*(1+s*s).^2))).^2;
    end
    f3(1,i)=k;
    i=i+1;
end
h=0.01;
x=0:h:5;
Y=diff(f3)/h;
% par le facon manuel
j=1;
for s=0:0.01:5
    k=0;
    for t=1:100
        f=a^2*exp(-((t-d).^2)/(2*(1+s*s).^2));
        k=k+(-4)*s/((1+s*s).^3)*((t-d).^2)*f*(sig_noisy(1,t)-f);
    end
    f3m(1,j)=k;
    j=j+1;
end
figure(13);

plot(x(:,1:length(Y)),Y,'c*',x,f3m,'k.');
legend('la fonction diff','la fonction manuelle');
title('par rapport ¨¤ ''s''');

% en utilisant les objets symboliques
syms a;
syms d;
syms s;
g = d;
for i=0:99   
    g=g+(sig_noisy(i+1)-a.^2.*exp(-(i-d).^2./(2.*(1+s.^2).^2))).^2;
end
g=g-d;
g(a,d,s)=g;
fa(a,d,s)=diff(g,a);
fd(a,d,s)=diff(g,d);
fs(a,d,s)=diff(g,s);
grad(a,d,s)=[fa(a,d,s) fd(a,d,s) fs(a,d,s)]'


%% Q15
xk=[2,4,15];%on initialise les valeurs et puis on calcul le module de ce 
% point pour entrer la boucle de la m¨¦thode des plus fortes pentes
[am,sm,dm,err]=deriveDir(xk(1,1),xk(1,2),xk(1,3));
k=0;
while(err>0.01)
    k=k+1;
    [dk(k,1),dk(k,2),dk(k,3),err]=deriveDir(xk(k,1),xk(k,2),xk(k,3));
    dk(k,:)=-dk(k,:);%la direction plus forte de decsente
    lambda=pas(dk(k,:),xk(k,:),0.01);
    %on stocke tous les it¨¦ration dans xk jusau'¨¤ l'erreur est acceptable
    xk(k+1,1)=xk(k,1)+lambda*dk(k,1);
    xk(k+1,2)=xk(k,2)+lambda*dk(k,2);
    xk(k+1,3)=xk(k,3)+lambda*dk(k,3);
end
points=[xk(k+1,:)]; 
point_final= -[xk(k+1,:)] %on s'arrete ¨¤ ce point.
valeur_FoctionCout=valeurFoctionCout(points(1,1),points(1,2),points(1,3));

xkdiff=diff(xk);
for i=1:size(xkdiff,1)
    mod_succ_x(i,1)=sqrt(xkdiff(i,1).^2+xkdiff(i,2).^2+xkdiff(i,3).^2);
end
% on calcul le module de d¨¦riv¨¦e de points calcul¨¦s pr¨¦c¨¦demment
figure(14);
plot(1:i,mod_succ_x(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦Èk-1||');title('l¡¯¨¦cart entre deux it¨¦r¨¦s successifs');

for i=1:size(xk,1)
    mod_optim(i,1)=sqrt((xk(i,1)-points(1,1)).^2+(xk(i,2)-points(1,2)).^2+(xk(i,3)-points(1,3)).^2);
end
figure(15);plot(1:i,mod_optim(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦È¡Ş ||');title('l¡¯¨¦cart ¨¤ l¡¯optimum ¦È¡Ş');

for i=1:size(xk,1)
    mod_fontion(i,1)=valeurFoctionCout(xk(i,1),xk(i,2),xk(i,3))-valeur_FoctionCout;
end
figure(16);plot(1:i,mod_fontion(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||f(¦Èk)-f(¦È¡Ş) ||');title(' l¡¯¨¦cart en terme de fonctions de co?t');

% la norme infinie du gradient ¡Î?¦ÈC (¦Èk)¡Î¡Ş
grad_points=[-dk(k,:)] %la norme infinie du gradient ¡Î?¦ÈC (¦Èk)¡Î¡Ş
% on s'arrete ¨¤ la position dont le gradient=(0.0048,0.0019,0.0069);

% la m¨¦thode d'objet symbolique
% [p_x p_y p_z]=fsxsteep(0.01,5,5,5);
% x=5;y=5;z=25;
% grad=grad(x,y,z);
% 
% faa(a,d,s)=diff(fa,a);
% fad(a,d,s)=diff(fa,d);
% fas(a,d,s)=diff(fa,s);
% 
% fda(a,d,s)=diff(fd,a);
% fdd(a,d,s)=diff(fd,d);
% fds(a,d,s)=diff(fd,s);
% 
% fsa(a,d,s)=diff(fs,a);
% fsd(a,d,s)=diff(fs,d);
% fss(a,d,s)=diff(fs,s);
% hessienne=([faa(x,y,z) fad(x,y,z) fas(x,y,z) ; fda(x,y,z) fdd(x,y,z) fds(x,y,z) ; fsa(x,y,z) fsd(x,y,z) fss(x,y,z)])';
% dir=-hessienne\grad(x,y,z);
% 
% x_k=[x y z]';
% e=0.01;
% while (norm(grad)>e)
%     pas=lambda(g,dir,x,y,z);
%     x_k=x_k+pas*dir; %ËÑË÷µ½µÄµã
%     x=[1 0 0]*x_k;
%     y=[0 1 0]*x_k;
%     z=[0 0 1]*x_k;
%     grad=[fa(x,y,z) fd(x,y,z) fs(x,y,z)]';
% end
%% Q16
% on part de diff¨¦rents points de d¨¦part.
xk=[2,1,25];%on initialise les valeurs et puis on calcul le module de ce 
% point pour entrer la boucle de la m¨¦thode des plus fortes pentes
[am,sm,dm,mod]=deriveDir(xk(1,1),xk(1,2),xk(1,3));
k=0;
while(mod>0.01)
    k=k+1;
    hess=matHessienne(xk(k,1),xk(k,2),xk(k,3),0.01);
    [L,tau]=factCholesky(hess);
    [d(k,1),d(k,2),d(k,3),mod]=deriveDir(xk(k,1),xk(k,2),xk(k,3));
    d(k,:)=-d(k,:);
    a=[d(k,1),d(k,2),d(k,3)];
    dk(k,:)=a/(L*L');
    lambda=pas(dk(k,:),xk(k,:),1);
    xk(k+1,1)=xk(k,1)+lambda*dk(k,1);
    xk(k+1,2)=xk(k,2)+lambda*dk(k,2);
    xk(k+1,3)=xk(k,3)+lambda*dk(k,3);
end
points=[xk(k+1,:)];
point_final= -[xk(k+1,:)] %on s'arrete ¨¤ ce point.
valeur_FoctionCout=valeurFoctionCout(points(1,1),points(1,2),points(1,3));

xkdiff=diff(xk);
for i=1:size(xkdiff,1)
    mod_succ_x(i,1)=sqrt(xkdiff(i,1).^2+xkdiff(i,2).^2+xkdiff(i,3).^2);
end
% on calcul le module de d¨¦riv¨¦e de points calcul¨¦s pr¨¦c¨¦demment
figure(17);
plot(1:i,mod_succ_x(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦Èk-1||');title('l¡¯¨¦cart entre deux it¨¦r¨¦s successifs');

for i=1:size(xk,1)
    mod_optim(i,1)=sqrt((xk(i,1)-points(1,1)).^2+(xk(i,2)-points(1,2)).^2+(xk(i,3)-points(1,3)).^2);
end
figure(18);plot(1:i,mod_optim(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦È¡Ş ||');title('l¡¯¨¦cart ¨¤ l¡¯optimum ¦È¡Ş');

for i=1:size(xk,1)
    mod_fontion(i,1)=valeurFoctionCout(xk(i,1),xk(i,2),xk(i,3))-valeur_FoctionCout;
end
figure(19);plot(1:i,mod_fontion(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||f(¦Èk)-f(¦È¡Ş) ||');title(' l¡¯¨¦cart en terme de fonctions de co?t');
grad_points=[-dk(k,:)]  %la norme infinie du gradient ¡Î?¦ÈC (¦Èk)¡Î¡Ş
% dans la question pr¨¦c¨¦demment,on part de point choisi apr¨¨s plusieurs
% essai de moi-meme qui vas nous pr¨¦senter le meilleur graphe d"it¨¦ration
% mais il n'a plus d'autre sens.
% maintenant,on part de point minimum calcul¨¦ pr¨¦c¨¦demment pour savoir le
% fonctionnement de points proches.
% En comparent les deux diff¨¦rents points choisis, on sais clairement que
% le fonctionnement est tr¨¨s diff¨¦rents par rapport aux points de d¨¦part.
%% Q17
xk=[2,4,15];%on initialise les valeurs et puis on calcul le module de ce 
% point pour entrer la boucle de la m¨¦thode de quasi-Newton 
[am,sm,dm,mod]=deriveDir(xk(1,1),xk(1,2),xk(1,3));
H0=inv(eye(3));
k=0;
while(mod>0.01)
    k=k+1;
    [d(k,1),d(k,2),d(k,3),mod1]=deriveDir(xk(k,1),xk(k,2),xk(k,3));
    dk(k,:)=-H0*(d(k,:)');
    lambda=pas(dk(k,:),xk(k,:),1);
    xk(k+1,1)=xk(k,1)+lambda*dk(k,1);
    xk(k+1,2)=xk(k,2)+lambda*dk(k,2);
    xk(k+1,3)=xk(k,3)+lambda*dk(k,3);
    [d(k+1,1),d(k+1,2),d(k+1,3),mod]=deriveDir(xk(k+1,1),xk(k+1,2),xk(k+1,3));
    yk=d(k+1,:)-d(k,:);
    difk=xk(k+1,:)-xk(k,:);
    H0=BFGS(yk',difk',H0);
end
points=[xk(k+1,:)];
point_final= -[xk(k+1,:)] %on s'arrete ¨¤ ce point.
valeur_FoctionCout=valeurFoctionCout(points(1,1),points(1,2),points(1,3));

xkdiff=diff(xk);
for i=1:size(xkdiff,1)
    mod_succ_x(i,1)=sqrt(xkdiff(i,1).^2+xkdiff(i,2).^2+xkdiff(i,3).^2);
end
% on calcul le module de d¨¦riv¨¦e de points calcul¨¦s pr¨¦c¨¦demment
figure(20);
plot(1:i,mod_succ_x(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦Èk-1||');title('l¡¯¨¦cart entre deux it¨¦r¨¦s successifs');

for i=1:size(xk,1)
    mod_optim(i,1)=sqrt((xk(i,1)-points(1,1)).^2+(xk(i,2)-points(1,2)).^2+(xk(i,3)-points(1,3)).^2);
end
figure(21);plot(1:i,mod_optim(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||¦Èk-¦È¡Ş ||');title('l¡¯¨¦cart ¨¤ l¡¯optimum ¦È¡Ş|');

for i=1:size(xk,1)
    mod_fontion(i,1)=valeurFoctionCout(xk(i,1),xk(i,2),xk(i,3))-valeur_FoctionCout;
end
figure(22);plot(1:i,mod_fontion(1:i),'k.');xlabel('nb d''i¨¦ration');ylabel('||f(¦Èk)-f(¦È¡Ş) ||');title(' l¡¯¨¦cart en terme de fonctions de co?t');
grad_points=[-dk(k,:)] %la norme infinie du gradient ¡Î?¦ÈC (¦Èk)¡Î¡Ş
% la m¨¦thode d'objet symbolique
% tao=0;
% hessienne=([faa(x,y,z) fad(x,y,z) fas(x,y,z) ; fda(x,y,z) fdd(x,y,z) fds(x,y,z) ; fsa(x,y,z) fsd(x,y,z) fss(x,y,z)]+tao*eye(3))';
% if(min(min(double(hessienne)))>0)
%     tao=0;
% else
%     tao=norm(double(hessienne));
% end
% hessienne=([faa(x,y,z) fad(x,y,z) fas(x,y,z) ; fda(x,y,z) fdd(x,y,z) fds(x,y,z) ; fsa(x,y,z) fsd(x,y,z) fss(x,y,z)]+tao*eye(3))';
%             
% while(norm(grad)<=0.01)  
%     d_k=-hessienne\grad;
%     grad_k_1=grad;
%     pas=lambda(g,d_k,x,y,z);
%     x_h=x_k+pas*d_k;
%     x=[1 0 0]*x_h;
%     y=[0 1 0]*x_h;
%     z=[0 0 1]*x_h;
%     
%     grad=[fa(x,y,z) fd(x,y,z) fs(x,y,z)]';
%     grad_k=grad;
%    
%     y_k_1 = grad_k-grad_k_1;
%     d_k_1 = x_h-x_k;
%     
%     hessienne = inv( (eye(3)-(d_k_1*y_k_1')/(d_k_1'*y_k_1)) *inv(hessienne)* (eye(3)-(y_k_1*d_k_1')/(d_k_1'*y_k_1)) + (d_k_1*d_k_1')/(d_k_1'*y_k_1) );
%     
% end

%% Q18
% on part ¨¤ diff¨¦rents points qui sont d¨¦j¨¤ utilis¨¦s dans la question n.16
% pour comparer diff¨¦rents points de d¨¦part. Alors on peut v¨¦rifier est que
% notre fonction marche corr¨¨ctemment.
x0=[2,4,15]
[x1,fval,flag]=fminunc('fonctionFminunc',x0)
x0=[2,1,25]
[x1,fval,flag]=fminunc('fonctionFminunc',x0)