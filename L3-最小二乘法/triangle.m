function[a]=triangle(A,b,n)
eps=0.01; %精度 
a=zeros(1,n); 
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
end
a(n)=b(n)/A(n,n);
for i=n-1:-1:1             %回代
    a(i)=b(i);
    for j=i+1:n
        a(i)=a(i)-A(i,j)*a(j);
    end
    a(i)=a(i)/A(i,i);
end



