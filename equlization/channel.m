function h=channel(L, P)
h=zeros(P,L);
for i=1:P
    if (i-1)/P>=0 && (i-1)/P<0.25
        h(i,1)=1;
    elseif (i-1)/P>=0.25 && (i-1)/P<0.5
        h(i,1)=-1;
    elseif (i-1)/P>=0.5 && (i-1)/P<0.75
        h(i,1)=1;
    else
        h(i,1)=-1;
    end
end
end



