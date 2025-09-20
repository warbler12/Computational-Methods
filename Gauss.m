n=2; %输入n
eps=0.01; %精度 
delta=1.0
A=[7 14.56;14.56 43.76]; 
b=[8.76 16.77];
x=zeros(1,n);            
for k=1:n-1             %消元_第k轮
    t=k;
    for i=k+1:n           %列主元
        if abs(A(i,k))>abs(A(t,k))
            t=i;
        end
    end
    if abs(A(t,k))<eps    %精度不符
        error("精度不符");
    end
    if(t~=k)            %交换
        for i=k:n
            a=A(k,i);A(k,i)=A(t,i);A(t,i)=a;
        end
        a=b(k);b(k)=b(t);b(t)=a;
    end
    for i=k+1:n
        A(i,k)=-A(i,k)/A(k,k);  %消元_第i行
        for j=k+1:n
            A(i,j)=A(i,j)+A(i,k)*A(k,j);
            
        end
        b(i)=b(i)+A(i,k)*b(k);
        A(i,k)=0;
    end
    delta=A(k,k)*delta;
end
x(n)=b(n)/A(n,n);
for i=n-1:-1:1                    %回代
    x(i)=b(i);
    for j=i+1:n
        x(i)=x(i)-A(i,j)*x(j);
    end
    x(i)=x(i)/A(i,i);
end
delta=A(n,n)*delta;
disp(A);
disp(b);
disp(x);
disp(delta);