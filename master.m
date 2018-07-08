clear all
clc
close all
%% signal model
M=5;%5 antennas
Delta=0.5;
N=20;
SNR=20;%in dB
theta=[-20 30];%two sources
f=[0.1 0.3];
[X,~,~]=gen_data(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
saveas(gcf,'initial.jpg');
% number of samples doubles
N=40;
[X,~,~]=gen_data(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
saveas(gcf,'sample_double.jpg');
% number of antennas doubles
N=20;
M=10;
[X,~,~]=gen_data(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
saveas(gcf,'antenna_double.jpg');
% angles between the sources becomes small
M=5;
theta=[0 5];
[X,~,~]=gen_data(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
saveas(gcf,'angle_small.jpg');
% frequency difference become small
theta=[-20 30];
f=[0 0.1];
[X,~,~]=gen_data(M,N,Delta,theta,f,SNR);
s=svd(X);
draw_svd(s,M,N,theta,f);
saveas(gcf,'frequency_difference_small.jpg');
%% estimation of directions
clear all
clc
close all
% MUSIC
M=5;%5 antennas
Delta=0.5;
N=20;
SNR=20;%in dB
theta=[-20 30];%two sources
f=[0.1 0.3];
d=2;
[X,A,S]=gen_data(M,N,Delta,theta,f,SNR);
theta_est=music(X,d,Delta);
saveas(gcf,'music.jpg');
% check the correctness in case there is no noise
theta_est_check=music(S,d,Delta);
saveas(gcf,'music_check.jpg');
%% estimation of frequencies
clear all
clc
close all
M=5;%5 antennas
Delta=0.5;
N=20;
SNR=20;%in dB
theta=[-20 30];%two sources
f=[0.1 0.3];
d=2;
[X,A,S]=gen_data(M,N,Delta,theta,f,SNR);
frequency_est=musicfreq(X,d);
saveas(gcf,'musicfreq.jpg');
% check the correctness in case there is no noise
frequency_est_check=musicfreq(S,d);
saveas(gcf,'musicfreq_check.jpg');
%% Comparison
