function f_est=musicfreq(X,d)
[~,~,V]=svd(X);
V=ctranspose(V);
[~,N]=size(X);
Vn=V(d+1:end,:);
f=0:0.01:1;
S=gen_f(f,N);
P=diag(S*ctranspose(S))./diag((S*ctranspose(Vn)*Vn*ctranspose(S)));%power in each direction
P=10*log10(abs(P));%rescale in dB
TF=islocalmax(P);
f_est=f(TF);
P_est=P(TF);
[~,I]=sort(P_est);
f_est=f_est(I(end-d+1:end));
P_est=P_est(I(end-d+1:end));
figure
plot(f, P, f_est, P_est, 'r*');
title('frequency estimation: MUSIC');
ylabel('Power[dB]');
xlabel('Normalized frequency[deg]');
end