function[output]=fs(A,i,k)
output=A(i,k);
    for j=1:k-1
        output=output-A(i,j)*A(j,k);
    end
end
