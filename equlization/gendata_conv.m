function [x,h]=gendata_conv(s,P,N,sigma,L)
s_ext=kron(s, [1; zeros(P-1,1)]);
h=channel(L,P);
x=conv(s_ext,h);
x=x(1:N*P);
v_N=sigma^2;
n=sqrt(v_N/2)*(randn(size(x))+1i*randn(size(x)));
x=x+n;
end