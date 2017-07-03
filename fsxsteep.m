function [x,y,z]=fsxsteep(e,x,y,z)
% x=fsxsteep(f,e,a,b)为输入函数 f为函数 e为允许误差 (a,b)为初始点;
global g;
syms a d s;
fa(a,d,s)=diff(g,a); %对x1求偏导数
fb(a,d,s)=diff(g,d); %对x2求偏导数
fc(a,d,s)=diff(g,s); %对x3求偏导数
grad=[fa(x,y,z) fb(x,y,z) fc(x,y,z)]'; %梯度值
x0=[x y z]';
while (norm(grad)>e)
    pas=lambda(g,-grad,x,y,z);
    x0=x0-grad.*pas; %搜索到的点
    x=[1 0 0]*x0;
    y=[0 1 0]*x0;
    z=[0 0 1]*x0;
    grad=[fa(x,y,z) fb(x,y,z) fc(x,y,z)]';
end