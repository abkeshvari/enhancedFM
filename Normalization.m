function [ D ] = Normalization( w,tau)
D=zeros(1,tau);
for i=1:tau
    if w(i)~=0
        D(i)=1/abs(w(i));
    else
        D(i)=1;
    end;
end;
end

