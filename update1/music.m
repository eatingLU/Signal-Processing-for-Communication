function theta_est=music(X,d,Delta)
[U, ~, ~]=svd(X);
[M,~]=size(X);
Un=U(:,d+1:end);
theta=-90:0.05:90;
A=gen_a(M,Delta,theta);
P=diag(ctranspose(A)*A)./diag((ctranspose(A)*Un*ctranspose(Un)*A));%power in each direction
P=10*log10(abs(P));%rescale in dB
TF=islocalmax(P);
theta_est=theta(TF);
P_est=P(TF);
if length(P_est)>=2
    [~,I]=sort(P_est);
    theta_est=theta_est(I(end-d+1:end));
    P_est=P_est(I(end-d+1:end));
    [theta_est,I]=sort(theta_est);
    P_est=P_est(I);
else
    theta_est=repmat(theta(TF),[2 1]);
    P_est=repmat(P(TF),[2 1]);
end
    
figure
plot(theta, P, theta_est, P_est, 'r*');
title('DOA estimation: MUSIC');
ylabel('Power[dB]');
xlabel('Angle[deg]');
str=[repmat('\theta:',2,1) num2str(theta_est')];
text(theta_est,P_est,cellstr(str))
end