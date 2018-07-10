function a=gen_a(M,Delta,theta)
gain=1;
d=length(theta);
theta_rad=theta/180*pi;% Convert degree into radian
a=zeros(M, d);
for i=1:M
    a(i, :)=exp(1j*2*pi*Delta*(i-1)*sin(theta_rad));
end
a=a.*gain;
end

