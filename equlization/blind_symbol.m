function s=blind_symbol(X,m,L)
% find the null-spaces
[~,~,V]=svd(X);
V_n=V(:,L+m:end);
V_n=[V_n; zeros(L+m-2,size(V_n,2))];% the first column of V_n,T
V_n_r=V_n;
for i=1:m+L-2
    V_n_r=[V_n_r circshift(V_n, i)];
end
[U_s, ~, ~]=svd(V_n_r);
s=U_s(:,end);%left null-sapce
end