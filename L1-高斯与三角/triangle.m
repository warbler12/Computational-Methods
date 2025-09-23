n=2; %输入n
eps=0.01; %精度 
delta=1.0;
A=[5 5327;5327 7277699]; %输入矩阵 
b=[271.4 369321.5
    
];
x=zeros(1,n);    
for k=1:n                 %第k轮
    s(k)=fs(A,k,k);       %L
    t=k;
    for i=k+1:n           %列主元
        s(i)=fs(A,i,k);
        if abs(s(i))>abs(s(t))
            t=i;
        end
    end
    if abs(s(t))<eps         %精度不符
        error("精度不符");
    end
    if(t~=k)              %交换
        a=s(t);s(t)=s(k);s(k)=a;
        for i=1:n    
            a=A(k,i);A(k,i)=A(t,i);A(t,i)=a;
        end
        a=b(k);b(k)=b(t);b(t)=a;
    end
    A(k,k)=s(k);
    for i=k+1:n
        A(i,k)=s(i)/s(k);
    end
    for i=k+1:n            %U
        A(k,i)=fu(A,k,i);
    end       
    for i=1:k-1             %y
        b(k)=b(k)-A(k,i)*b(i);
    end
    delta=abs(A(k,k))*delta;
end
x(n)=b(n)/A(n,n);
for i=n-1:-1:1             %回代
    x(i)=b(i);
    for j=i+1:n
        x(i)=x(i)-A(i,j)*x(j);
    end
    x(i)=x(i)/abs(A(i,i));
end
disp(A);
disp(b);
disp(x);
disp(delta);

function[output]=fs(A,i,k)
output=A(i,k);
    for j=1:k-1
        output=output-A(i,j)*A(j,k);
    end
end
function[output]=fu(A,k,j)
output=A(k,j);
    for i=1:k-1
        output=output-A(k,i)*A(i,j);
    end
end
