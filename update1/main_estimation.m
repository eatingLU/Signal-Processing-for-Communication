clear all
clc
close all
% signal model
M=5;%5 antennas
Delta=0.5;
N=20;
SNR=20;%in dB
theta=[-20 30];%two sources
f=[0.1 0.3];
[X,~,~]=gendata(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
% saveas(gcf,'initial.jpg');
% number of samples doubles
N=40;
[X,~,~]=gendata(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
% saveas(gcf,'sample_double.jpg');
% number of antennas doubles
N=20;
M=10;
[X,~,~]=gendata(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
% saveas(gcf,'antenna_double.jpg');
% angles between the sources becomes small
M=5;
theta=[0 5];
[X,~,~]=gendata(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
% saveas(gcf,'angle_small.jpg');
% frequency difference become small
theta=[-20 30];
f=[0 0.1];
[X,~,X_clean,S]=gendata(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
% saveas(gcf,'frequency_difference_small.jpg');
%% estimation of directions
% MUSIC
close all 
clc
clear all
M=5;%5 antennas
Delta=0.5;
N=20;
SNR=20;%in dB
theta=[-20 30];%two sources
f=[0.1 0.3];
d=2;
[X,~,X_clean,~]=gendata(M,N,Delta,theta,f,SNR);
theta_est=music(X,d,Delta);
% saveas(gcf,'music.jpg');
% check the correctness in case there is no noise
theta_est_check=music(X_clean,d,Delta);
% saveas(gcf,'music_check.jpg');
%% estimation of frequencies
M=5;
Delta=0.5;
N=20;
SNR=20;
theta=[-20 30];
f=[0.1 0.3];
d=2;
[X,~,X_clean,~]=gendata(M,N,Delta,theta,f,SNR);
frequency_est=musicfreq(X,d);
% saveas(gcf,'musicfreq.jpg');
% check the correctness in case there is no noise
frequency_est_check=musicfreq(X_clean,d);
% (gcf,'musicfreq_check.jpg');
%% estimation performance
SNR=0:4:20;
d=2;
M=3;
N=20;
theta=[-20 30];
f=[0.1 0.12];
Delta=0.5;
rep=1000;
theta_est=zeros(2,rep);
frequency_est=zeros(2,rep);
std_theta=zeros(2,length(SNR));
std_frequency=zeros(2,length(SNR));
m_theta=zeros(2,length(SNR));
m_frequency=zeros(2,length(SNR));
for i=1:length(SNR)
    for j=1:rep
        [X,A,S]=gendata(M,N,Delta,theta,f,SNR(i));
        theta_est(:,j)=music(X,d,Delta);
        frequency_est(:,j)=musicfreq(X,d);
    end
    m_theta(:,i)=mean(theta_est,2);
    m_frequency(:,i)=mean(frequency_est,2);
    std_theta(:,i)=sqrt(mean((theta_est-repmat(theta',[1 rep])).^2,2));
    std_frequency(:,i)=sqrt(mean((frequency_est-repmat(f',[1 rep])).^2,2));
end
figure
subplot(221)
plot(SNR,m_theta(1,:),SNR,m_theta(2,:),SNR,repmat(theta(1),[1 length(SNR)]),...
    'b--',SNR,repmat(theta(2),[1 length(SNR)]),'r--');
title('estimation performance of DOA (mean values)');
xlabel('SNR[dB]');
ylabel('DOA[deg]');
legend('estimated source 1','estimated source 2','source 1','source 2');
subplot(222)
plot(SNR,m_frequency(1,:),SNR,m_frequency(2,:),SNR,repmat(f(1),[1 length(SNR)]),...
    'b--',SNR,repmat(f(2),[1 length(SNR)]),'r--');
title('estimation performance of frequencies (mean values)');
xlabel('SNR[dB]');
ylabel('normalized frequency');
legend('estimated source 1','estimated source 2','source 1','source 2');
subplot(223)
plot(SNR,std_theta(1,:),SNR,std_theta(2,:));
title('estimation performance of DOA (stand deviations)');
xlabel('SNR[dB]');
ylabel('stand deviations');
legend('source 1','source 2');
subplot(224)
plot(SNR,std_frequency(1,:),SNR,std_frequency(2,:));
title('estimation performance of frequencies (stand deviations)');
xlabel('SNR[dB]');
ylabel('stand deviations');
legend('source 1','source 2');
% saveas(gcf,'performance.jpg');
%% zero-forcing beamformer
% based on direction estimates from music
% with A known
M=5;
Delta=0.5;
N=20;
SNR=20;
theta=[-20 30];
f=[0.1 0.12];
d=2;
[X,~,~,S]=gendata(M,N,Delta,theta,f,SNR);
theta_est=music(X,d,Delta);
% saveas(gcf,'theta_est.jpg');
A_known=gen_a(M, Delta, theta_est);
W=pinv(A_known);
S_est=W*X;
figure
subplot(211)
plot(real(S_est(1,:)),imag(S_est(1,:)),'bx',real(S(1,:)),imag(S(1,:)),'ro');
title('ZF beamformer based on direction estimates for source 1');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 1','source 1');
axis([-1.5 1.5 -1.5 1.5]);
subplot(212)
plot(real(S_est(2,:)),imag(S_est(2,:)),'bx',real(S(2,:)),imag(S(2,:)),'ro');
title('ZF beamformer based on direction estimates for source 2');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 2','source 2');
axis([-1.5 1.5 -1.5 1.5]);
% saveas(gcf,'direction_based_ZF.jpg');
% based on frequency estimate from 'musicfreq'
% with S known
f_est=musicfreq(X,d);
saveas(gcf,'f_est.jpg');
S_known=gen_f(f_est,N);
W=pinv(X*pinv(S_known));
S_est=W*X;
figure
subplot(211)
plot(real(S_est(1,:)),imag(S_est(1,:)),'bx',real(S(1,:)),imag(S(1,:)),'ro');
title('ZF beamformer based on frequency estimates for source 1');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 1','source 1');
axis([-1.5 1.5 -1.5 1.5]);
subplot(212)
plot(real(S_est(2,:)),imag(S_est(2,:)),'bx',real(S(2,:)),imag(S(2,:)),'ro');
title('ZF beamformer based on frequency estimates for source 2');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 2','source 2');
axis([-1.5 1.5 -1.5 1.5]);
% saveas(gcf,'frequency_based_ZF.jpg');
%% test the correctness
M=5;
Delta=0.5;
N=20;
SNR=20;
theta=[-20 30];
f=[0.1 0.12];
d=2;
[~,~,X_clean,S]=gendata(M,N,Delta,theta,f,SNR);
theta_est=music(X_clean,d,Delta);
saveas(gcf,'correct_theta_est.jpg');
A_known=gen_a(M, Delta, theta_est);
W=pinv(A_known);
S_est=W*X_clean;
figure
subplot(211)
plot(real(S_est(1,:)),imag(S_est(1,:)),'bx',real(S(1,:)),imag(S(1,:)),'ro');
title('correctness test (ZF beamformer based on direction estimates for source 1)');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 1','source 1');
axis([-1.5 1.5 -1.5 1.5]);
subplot(212)
plot(real(S_est(2,:)),imag(S_est(2,:)),'bx',real(S(2,:)),imag(S(2,:)),'ro');
title('correctness test (ZF beamformer based on direction estimates for source 2)');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 2','source 2');
axis([-1.5 1.5 -1.5 1.5]);
% saveas(gcf,'correctness_direction_based_ZF.jpg');
% based on frequency estimate from 'musicfreq'
% with S known
f_est=musicfreq(X_clean,d);
saveas(gcf,'correct_f_est.jpg');
S_known=gen_f(f_est,N);
W=pinv(X_clean*pinv(S_known));
S_est=W*X_clean;
figure
subplot(211)
plot(real(S_est(1,:)),imag(S_est(1,:)),'bx',real(S(1,:)),imag(S(1,:)),'ro');
title('correctness test (ZF beamformer based on frequency estimates for source 1)');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 1','source 1');
axis([-1.5 1.5 -1.5 1.5]);
subplot(212)
plot(real(S_est(2,:)),imag(S_est(2,:)),'bx',real(S(2,:)),imag(S(2,:)),'ro');
title('correctness test (ZF beamformer based on frequency estimates for source 2)');
xlabel('real(S_{est})');
ylabel('imaginary(S_{est})');
legend('estimated source 2','source 2');
axis([-1.5 1.5 -1.5 1.5]);
% saveas(gcf,'correctness_frequency_based_ZF.jpg');
%% spatial response
M=5;
Delta=0.5;
N=20;
SNR=10;
theta=[-20 30];
f=[0.1 0.12];
d=2;
[X,~,~,S]=gendata(M,N,Delta,theta,f,SNR);
% ZF beamformer based on direction estimates
theta_est=music(X,d,Delta);
% saveas(gcf,'theta_est_SNR10.jpg');
A_known=gen_a(M, Delta, theta_est);
W=pinv(A_known);
theta=-90:90;
a=gen_a(M,Delta,theta);
y=abs(W*a);
figure
plot(theta,y(1,:),theta,y(2,:));
title('zero-forcing beamformer spatial response(based on direction estimates)');
xlabel('angle[deg]');
legend('source 1','source 2');
% saveas(gcf,'spatial_response_direction_based_ZF.jpg');
% ZF beamformer based on frequency estimates
f_est=musicfreq(X,d);
% saveas(gcf,'f_est_SNR10.jpg');
S_known=gen_f(f_est,N);
W=pinv(X*pinv(S_known));
theta=-90:90;
a=gen_a(M,Delta,theta);
y=abs(W*a);
figure
plot(theta,y(1,:),'-.o',theta,y(2,:),'-.x');
title('zero-forcing beamformer spatial response(based on frequency estimates)');
xlabel('angle[deg]');
legend('source 1','source 2');
% saveas(gcf,'spatial_response_frequency_based_ZF.jpg');
%% CMA algorithm
M=5;
Delta=0.5;
N=5000;
SNR=20;
theta=[-20 30];
f=[0.1 0.3];
d=2;
[X,~,X_clean,~]=gendata(M,N,Delta,theta,f,SNR);
w_init=ones(M,1);
mu=0.001;
[w,m]=cma(X,mu,w_init,N);
y=ctranspose(w)*X;
c_x=0;c_y=0;r=1;
figure
subplot(121)
plot(X(1,:),'.');
hold on
rectangle('Position',[c_x-r,c_y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1);
hold off
ylabel('imag(y)');
xlabel('real(y)');
title('before beamforming');
axis equal;
subplot(122)
plot(y(1,:),'.');
hold on
rectangle('Position',[c_x-r,c_y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1);
hold off
ylabel('imag(y)');
xlabel('real(y)');
title('after beamforming');
axis equal
% saveas(gcf,'CMA_beamforming_result.jpg');
figure
plot(m);
title('modulus of output');
xlabel('time(samples)');
ylabel('abs(y)');
% saveas(gcf,'CMA_convergence.jpg');
% varify the correctness of the algorithm
[w,m]=cma(X_clean,mu,w_init,N);
y=ctranspose(w)*X_clean;
theta=-90:90;
a=gen_a(M,Delta,theta);
y_sr=abs(ctranspose(w)*a);
c_x=0;c_y=0;r=1;
figure
subplot(121)
plot(X_clean(1,:),'*');
hold on
rectangle('Position',[c_x-r,c_y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1);
hold off
ylabel('imag(y)');
xlabel('real(y)');
title('before beamforming');
axis equal;
subplot(122)
plot(y(1,:),'*');
hold on
rectangle('Position',[c_x-r,c_y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1);
hold off
ylabel('imag(y)');
xlabel('real(y)');
title('after beamforming');
axis equal
% saveas(gcf,'varify_CMA_beamforming_result.jpg');
figure
plot(m);
title('modulus of output');
xlabel('time(samples)');
ylabel('abs(y)');
% saveas(gcf,'varify_CMA_convergence.jpg');
figure
plot(theta,y_sr);
title('CMA beamfomer spatial response');
xlabel('angle[deg]');
% saveas(gcf,'CMA_spatial_response.jpg');
%% random initial weight vectors
close all
clc
clear all
M=5;
Delta=0.5;
N=5000;
SNR=20;
theta=[-20 30];
f=[0.1 0.3];
d=2;
[X,~,X_clean,~]=gendata(M,N,Delta,theta,f,SNR);
w_init=[1 0 0 0; 1 1 0 0; 1 0 1 0; 1 0 0 1; 1 0 0 0];
% w_init=eye(M);
mu=0.001;
theta=-90:90;
a=gen_a(M,Delta,theta);
for i=1:4
    [w,m]=cma(X,mu,w_init(:,i),N);
    y_sr=abs(ctranspose(w)*a);
    subplot(2,2,i)
    plot(theta,y_sr);
    title(['CMA beamfomer spatial response with w_{init}=[' num2str(w_init(:,i)') ']']);
    xlabel('angle[deg]');
end
%saveas(gcf,'CMA_spatial_response_random.jpg');
