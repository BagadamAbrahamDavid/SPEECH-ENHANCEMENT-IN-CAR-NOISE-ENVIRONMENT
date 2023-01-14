function y=rfft(x,n,d)
% 0002 %RFFT     Calculate the DFT of real data Y=(X,N,D)
% 0003 % Data is truncated/padded to length N if specified.
% 0004 %   N even:    (N+2)/2 points are returned with
% 0005 %             the first and last being real
% 0006 %   N odd:    (N+1)/2 points are returned with the
% 0007 %             first being real
% 0008 % In all cases fix(1+N/2) points are returned
% 0009 % D is the dimension along which to do the DFT
% 0010 
% 0011 
% 0012 

 s=size(x);
 if prod(s)==1
     y=x
 else
     if nargin <3 || isempty(d)
         d=find(s>1,1);
         if nargin<2
             n=s(d);
         end
     end
     if isempty(n) 
         n=s(d);
     end
    y=fft(x,n,d);
     y=reshape(y,prod(s(1:d-1)),n,prod(s(d+1:end))); 
     s(d)=1+fix(n/2);
     y(:,s(d)+1:end,:)=[];
     y=reshape(y,s);
 end