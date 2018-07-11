function X=construct_X(x,P,N,m)
X=zeros(m*P,N-m+1);
for i=1:m*P
    for j=1:N-m+1
        X(i,j)=x(i+P*(j-1));
    end
end
end