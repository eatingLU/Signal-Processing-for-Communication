function s=source(N)
QPSKmap=[1, 1i, -1, -1i];        
temp=randi(4,N,1);% generate N by 1 random integers in interval[1,4]
s=QPSKmap(temp)';
end