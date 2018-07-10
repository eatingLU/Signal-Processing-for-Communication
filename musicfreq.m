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
[f_est,I]=sort(f_est);
P_est=P_est(I);
if f_est(1)<0.2*f_est(2)
    f_est=repmat(f_est(2),[1 2]);
    P_est=repmat(P_est(2),[1 2]);
    figure
    plot(f, P, f_est, P_est, 'r*');
    title('frequency estimation: MUSIC');
    ylabel('Power[dB]');
    xlabel('Normalized frequency[deg]');
    str=[repmat('f:',2,1) num2str(f_est')];
    text(f_est,P_est,cellstr(str))
else
    figure
    plot(f, P, f_est, P_est, 'r*');
    title('frequency estimation: MUSIC');
    ylabel('Power[dB]');
    xlabel('Normalized frequency[deg]');
    str=[repmat('f:',2,1) num2str(f_est')];
    text(f_est,P_est-[0 P_est(2)/20]',cellstr(str))
end
end