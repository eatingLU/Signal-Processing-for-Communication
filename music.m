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
[~,I]=sort(P_est);
theta_est=theta_est(I(end-d+1:end));
P_est=P_est(I(end-d+1:end));
figure
plot(theta, P, theta_est, P_est, 'r*');
title('DOA estimation: MUSIC');
ylabel('Power[dB]');
xlabel('Angle[deg]');
end