function []=draw_svd(s,M,N,theta,f)
figure
stem(s);
ylabel('singular value');
xlabel('index');
title(['SVD of data matrix with M=' num2str(M) ',N=' num2str(N) ',\theta=[' ...
num2str(theta(1)) ',' num2str(theta(2)) '] and f=[' num2str(f(1)) ',' num2str(f(2)) ']']);
end