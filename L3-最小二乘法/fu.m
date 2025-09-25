function[output]=fu(A,k,j)
output=A(k,j);
    for i=1:k-1
        output=output-A(k,i)*A(i,j);
    end
end