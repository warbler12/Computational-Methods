n=10;                    %输入层数n
eps=0.5*0.000001;              %精度
x=zeros(n,2^n+1);
x(1,1)=0.1^100;x(1,2)=1;      %输入积分区间
h=x(1,2)-x(1,1);            %初步长
T=zeros(n,4);    
for k=1:n                 %第k层
    if k==1                                 %T
        T(1,1)=(f(x(1,1))+f(x(1,2)))*h/2;
    else    
        T(k,1)=0.5*T(k-1,1);
        for i=2:2:2^(k-1)+1
            T(k,1)=T(k,1)+h*f(x(k,i));
        end
        m=min(k,4);                             %S/C/R
        for i=2:m
            T(k,i)=(4^(i-1)*T(k,i-1)-T(k-1,i-1))/(4^(i-1)-1);
        end
        if m==k                                 %delta
            delta=T(k,m)-T(k-1,m-1);
        else
            delta=T(k,m)-T(k-1,m);
        end 
        if abs(delta)<eps 
            break;
        end
    end
    if k<n
        x(k+1,1)=x(1,1);
        for i=2:2^(k-1)+1                       %x
            x(k+1,2*i-1)=x(k,i);
            x(k+1,2*(i-1))=(x(k,i-1)+x(k,i))/2;
        end
    end
    h=h/2;
end
if abs(delta)>eps
    error("精度不符");
else
    fprintf("T=\n");
    for i=1:k
        fprintf("\n");
        for j=1:4
            fprintf('   %5.7f  ',T(i,j));
        end
    end
    delta=delta
end

function output =f(x)   %输入被积函数
    output=sin(x)/x;
end

