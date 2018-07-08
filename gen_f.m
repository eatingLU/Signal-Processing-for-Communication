function S=gen_f(f,N)
d=length(f);
S=zeros(d,N);
for i=1:d
    for j=1:N
        S(i,j)=exp(1i*2*pi*f(i)*j);
    end
end
end