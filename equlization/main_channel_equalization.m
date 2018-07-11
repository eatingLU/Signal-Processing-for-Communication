close all
clc
clear all
P=4;
N=500;
sigma=0;
L=1;
s=source(N+L-1);
[x,~]=gendata_conv(s,P,N,sigma,L);
X=construct_X(x,P,N);
rX=rank(X);
% double P
P=8;
x=gendata_conv(s,P,N,sigma,L);
X=construct_X(x,P,N);
rX_double=rank(X);
%% symbol estimation
% zero-forcing equalizer
P=4;
N=500;
sigma=0.5;
L=1;
s=source(N+L-1);
[x,h]=gendata_conv(s,P,N,sigma,L);
X=construct_X(x,P,N);
H=[h,zeros(size(h));zeros(size(h)),h];
e=[1;0];
w=pinv(H*ctranspose(H))*H*e;
s_est=ctranspose(w)*X;
figure
subplot(221)
plot(s_est,'.');
title('estimated source symbols using ZF beamformer');
xlabel('real(s_{est})');
ylabel('imaginary(s_{est})');
axis equal
% wiener equalizer
w=inv(H*ctranspose(H)+sigma^2*eye(2*P))*H*e;
s_est=ctranspose(w)*X;
subplot(223)
plot(s_est,'.');
title('estimated source symbols using Wiener beamformer');
xlabel('real(s_{est})');
ylabel('imaginary(s_{est})');
axis equal
% double P
P=8;
[x,h]=gendata_conv(s,P,N,sigma,L);
X=construct_X(x,P,N);
H=[h,zeros(size(h));zeros(size(h)),h];
e=[1;0];
w=pinv(H*ctranspose(H))*H*e;
s_est=ctranspose(w)*X;
subplot(222)
plot(s_est,'.');
title('estimated source symbols using ZF beamformer (doubled P)');
xlabel('real(s_{est})');
ylabel('imaginary(s_{est})');
axis equal
% wiener equalizer
w=inv(H*ctranspose(H)+sigma^2*eye(size(H*ctranspose(H))))*H*e;
s_est=ctranspose(w)*X;
subplot(224)
plot(s_est,'.');
title('estimated source symbols using Wiener beamformer (doubled P)');
xlabel('real(s_{est})');
ylabel('imaginary(s_{est})');
axis equal
saveas(gcf,'symbol_estimation.jpg');
%% channel estimation
P=4;
N=500;
sigma=0.5;
L=1;
s=source(N+L-1);
[x,h]=gendata_conv(s,P,N,sigma,L);
index=0:1/P:L-1/P;
figure
subplot(421)
stem(index,real(h));
title('real part of channel(ground truth)');
ylabel('real(h)');
xlabel('index');
axis([0 L -1.2 1.2]);
subplot(422)
stem(index,imag(h));
title('imaginary part of channel(ground truth)');
ylabel('imaginary(h)');
xlabel('index');
axis([0 L -1.2 1.2]);
%channel estimation
L=[1,2,3];
for i=1:length(L)
    h_hat=channel_estimator(x,s,L(i),P);
    index=0:1/P:L(i)-1/P;
    subplot(4,2,i*2+1)
    stem(index,real(h_hat));
    title(['real part of channel(L=' num2str(L(i)) ')']);
    ylabel('real(h)');
    xlabel('index');
    axis([0 L(i) -1.2 1.2]);
    subplot(4,2,i*2+2)
    stem(index,imag(h_hat));
    title(['imaginary part of channel(L=' num2str(L(i)) ')']);
    ylabel('imaginary(h)');
    xlabel('index');
    axis([0 L(i) -1.2 1.2]);
end
% saveas(gcf,'channel_estimation.jpg');
%% blind spatial processing
P=4;
N=500;
sigma=0.5;
L=1;
s=source(N+L-1);
[x,h]=gendata_conv(s,P,N,sigma,L);
m=2;% 2 symbol periods
X=construct_X(x,P,N,m);
% compute SVD
[U,S,V]=svd(X);
Un=U(:,m+L:end);% rank of S is m+L, X=HS
Vn=V(:,m+L:end);
%% blind channel estimation
P=4;
N=500;
sigma=0.5;
L=1;
s=source(N+L-1);
[x,h]=gendata_conv(s,P,N,sigma,L);
m=2;% 2 symbol periods
X=construct_X(x,P,N,m);
h_hat=blind_channel(X,m,L);
index=0:1/P:L-1/P;
figure
% true channel
subplot(221)
stem(index,real(h));
title('real part of channel(ground truth)');
ylabel('real(h)');
xlabel('index');
axis([0 L -1.2 1.2]);
subplot(222)
stem(index,imag(h));
title('imaginary part of channel(ground truth)');
ylabel('imaginary(h)');
xlabel('index');
axis([0 L -1.2 1.2]);
% estimated channel
subplot(223)
stem(index,real(h_hat));
title('real part of channel(estimated)');
ylabel('real(h_{hat})');
xlabel('index');
axis([0 L -1.2 1.2]);
subplot(224)
stem(index,imag(h_hat));
title('imaginary part of channel(estimated)');
ylabel('imaginary(h_{hat})');
xlabel('index');
axis([0 L -1.2 1.2]);
% saveas(gcf,'blind_channel.jpg');
%% estimated symbol
close all
clc
clear all
P=4;
N=500;
sigma=0.5;
L=1;
s=source(N+L-1);
[x,~]=gendata_conv(s,P,N,sigma,L);
m=2;% 2 symbol periods
X=construct_X(x,P,N,m);
s_hat=blind_symbol(X,m,L);
index=0:1/P:L-1/P;
figure
% true symbol
subplot(121)
plot(s,'*');
title('true symbol');
xlabel('real(s)');
ylabel('imaginary(s)');
axis([-1.1 1.1 -1.1 1.1]);
% estimated symbol
subplot(122)
plot(s_hat,'.');
title('estimated symbol');
xlabel('real(s_{hat})');
ylabel('imaginary(s_{hat})');
axis equal
saveas(gcf,'blind_symbol.jpg');