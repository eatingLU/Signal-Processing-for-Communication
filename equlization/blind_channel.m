function h=blind_channel(X,m,L)
P=size(X,1)/m;
% find the null-spaces
[U,~,~]=svd(X);
U_n=U(:,L+m:end)';
% find the right null-space of U_n
% reconstruct U_n matrix in size (m+L-1)*P
for i=1:m
    if i==1
        U_n_r=U_n(:,(i-1)*P+1:i*P);
    elseif i==m
        U_n_r=[U_n_r; U_n(:,(i-1)*P+1:i*P); zeros(size(U_n,1)*(L-1),P)];
    else
        U_n_r=[U_n_r; U_n(:,(i-1)*P+1:i*P)];
    end
end
%reconstruct U_n matrix in size (m+L-1)*P*L
for i=1:L-1
    U_n_r=[U_n_r circshift(U_n_r, size(U_n,1))];
end
%perfoem SVD to find the null-space of U_n
[~,~,V_h]=svd(U_n_r);
h=V_h(:,end);%right null-space
end