x1=[2 4 6 8]; 
y=[2 11 28 40];
plot(x1,y);
hold on
%一次
n=2;                %φ个数（次数）+1（即矩阵行数）
nx=4;               %x个数
A=zeros(n,n);       %(φ内积矩阵)
b=zeros(1,n);       %等号右边含y列向量
x=zeros(n*2-1,nx);     
y0=zeros(1,nx);
s=0;
for i=1:n*2-1        %计算x的幂
    for j=1:nx
        x(i,j)=x1(j)^(i-1);
    end
end
for i=1:n           %计算A,b
    for j=i:n
        for m=1:nx
            A(i,j)=A(i,j)+x(i+j-1,m);
            if j==i
                b(i)=b(i)+y(m)*x(i,m);
            end
        end
    end
    A(j,i)=A(i,j);
end
a=triangle(A,b,n);  %ak
for i=1:nx          %拟合得到
    y0(i)=x1(i)*a(2)+a(1);
end
delta=y-y0;     %最大偏差
for i=1:nx
    delta(i)=abs(delta(i));
end
deltamax=max(delta);
for i=1:nx
    s=s+delta(i)^2;
end
s=s^0.5;
fprintf("deltamax=%f\n",deltamax);
fprintf("s=%f\n",s);
fprintf("y=(%f)*x+(%f)\n",a(2),a(1));
plot(x1,y0);
hold off