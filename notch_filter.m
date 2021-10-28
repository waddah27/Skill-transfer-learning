function b=notch_filter(fc,T1,Fs)
 
N=floor(T1*Fs/2); 
t=[-N:N]/Fs; 
h=zeros(size(t)); 
for i=1:length(fc) 
h0=cos(2*pi*t*fc(i)); 
h1=h0.*(1+cos(pi/N*[-N:N]))/2; 
h1=h1/sum(h1.*h0); 
h=h+h1;
end 
b=zeros(size(h)); 
b(N+1)=1; 
b=b-h; 
end

