function [X,A,S]=gen_data(M,N,Delta,theta,f,SNR)
% for simplicity, assume signal variance is equal to 1
% compute array response vector
d=length(theta);%number of sources
a=zeros(M,d);
r=theta*pi/180;%convert to radians
for i=1:M
    for j=1:d
        a(i,j)=exp(1i*2*pi*(i-1)*Delta*sin(r(j)));
    end
end
% source signal
S=zeros(d,N);
for i=1:d
    for j=1:N
        S(i,j)=exp(1i*2*pi*f(i)*j);
    end
end
v_S=[1 1];
% compute the noise variance
v_N=v_S/10^(SNR/10);
% compute the noise
X=zeros(M,N,d);
X_clean=zeros(M,N,d);
for i=1:d
    n=sqrt(v_N(i)/2)*(randn(M,N)+1i*randn(M,N));
    X_clean(:,:,i)=a(:,i)*S(i,:);
    X(:,:,i)=X_clean(:,:,i)+n;
end
X=sum(X,3);
S=sum(X_clean,3);
A=a;
end
