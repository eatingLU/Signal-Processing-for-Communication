function [w,m]=cma(X,mu,w_init,N)
w=w_init;
m=zeros(N,1);
for i=1:N
    k=max(1,min(5000,floor(5000*rand)));
    y=ctranspose(w)*X(:,k);
    m(i)=abs(y);
    w=w-mu*X(:,k)*(abs(y)^2-1)*conj(y);
end
end