function h=channel_estimator(x,s,L,P)
N=length(s)-L+1;
% construct S matrix
i=length(s);
S=zeros(N, L);
j=1;
while(i>=L)
    S(j,:)=s(i:-1:i-L+1);
    i=i-1;
    j=j+1;
end
% construct x vector
x=x(1:N*P);
x_r=flip(x);
for i=1:P:N*P
    x_r(i:i+P-1)=flip(x_r(i:i+P-1));
end
% pilot estimation
I=eye(P);
h=(kron(pinv(S),I))*x_r;
end