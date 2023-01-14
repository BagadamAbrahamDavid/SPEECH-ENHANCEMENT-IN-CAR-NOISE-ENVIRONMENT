 function x=irfft(y,n,d)
 s=size(y);
 ps=prod(s);
 ns=length(s);
 if ps==1
     x=y;
 else
     if nargin <3 || isempty(d)
         d=find(s>1,1);
     end
     m=s(d);
     k=ps/m;     % number of fft's to do
     if d==1
         v=reshape(y,m,k);
     else
         v=reshape(permute(y,[d:ns 1:d-1]),m,k);
     end
     if nargin<2 || isempty(n)
         n=2*m-2;        % default output length
     else
         mm=1+fix(n/2);          % expected input length
         if mm>m v=[v; zeros(mm-m,k)];   % zero pad
      elseif mm<m v(mm+1:m,:)=[];     % or truncate
         end
        m=mm;
    end
  if rem(n,2)        % odd output length
      x=real(ifft([v;conj(v(m:-1:2,:))],[],1));    % do it the long way
  else            % even output length
         v(m,:)=real(v(m,:));    % force nyquist element real
      w=ones(1,k);
      %  t=[cumprod([-0.5i; exp(2i*pi/n)*ones(m-2,1)]); 0.5i];
      t=-0.5i* exp((2i*pi/n)*(0:m-1)).';
         z=(t(:,w)+0.5).*(conj(flipud(v))-v)+v;
         z(m,:)=[];
         zz=ifft(z,[],1);
        x=zeros(n,k);
         x(1:2:n,:)=real(zz);
         x(2:2:n,:)=imag(zz);
     end
  s(d)=n;         % change output dimension
     if d==1
         x=reshape(x,s);
     else
         x=permute(reshape(x,s([d:ns 1:d-1])),[ns+2-d:ns 1:ns+1-d]);
     end
 end