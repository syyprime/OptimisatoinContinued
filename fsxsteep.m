function [x,y,z]=fsxsteep(e,x,y,z)
% x=fsxsteep(f,e,a,b)Ϊ���뺯�� fΪ���� eΪ������� (a,b)Ϊ��ʼ��;
global g;
syms a d s;
fa(a,d,s)=diff(g,a); %��x1��ƫ����
fb(a,d,s)=diff(g,d); %��x2��ƫ����
fc(a,d,s)=diff(g,s); %��x3��ƫ����
grad=[fa(x,y,z) fb(x,y,z) fc(x,y,z)]'; %�ݶ�ֵ
x0=[x y z]';
while (norm(grad)>e)
    pas=lambda(g,-grad,x,y,z);
    x0=x0-grad.*pas; %�������ĵ�
    x=[1 0 0]*x0;
    y=[0 1 0]*x0;
    z=[0 0 1]*x0;
    grad=[fa(x,y,z) fb(x,y,z) fc(x,y,z)]';
end