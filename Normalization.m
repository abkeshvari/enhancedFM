function [ D ] = Normalization( w,tau)
%%% Please cite the paper properly if you use the code. 
%%% "Keshvari, Abolfazl. 2016. An Enhanced Fourier-Motzkin Method for DEA."
D=zeros(1,tau);
for i=1:tau
    if w(i)~=0
        D(i)=1/abs(w(i));
    else
        D(i)=1;
    end;
end;
end

