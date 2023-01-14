function [f,t,w]=enframe(x,win,inc,m)
% 0002 %ENFRAME split signal up into (overlapping) frames: one per row. [F,T]=(X,WIN,INC)
% 0003 %
% 0004 % Usage:  (1) f=enframe(x,n)     % split into frames of length n
% 0005 %         (2) f=enframe(x,hamming(n,'periodic'),n/4)     % use a 75% overlapped Hamming window of length n
% 0006 %         (3) calculate spectrogram in units of power per Hz
% 0007 %
% 0008 %               W=hamming(NW);                      % analysis window (NW = fft length)
% 0009 %               W=W/sqrt(FS*sum(W.^2));             % normalize to give power per Hz (FS = sample freq)
% 0010 %               P=rfft(enframe(S,W,INC);,nfft,2);   % computer first half of fft (INC = frame increment in samples)
% 0011 %               P(:,2:end-1)=2*P(:,2:end-1);        % double to account for -ve frequencies (except DC and Nyquist)
% 0012 %
% 0013 %         (3) frequency domain frame-based processing:
% 0014 %
% 0015 %               S=...;                              % input signal
% 0016 %               OV=2;                               % overlap factor of 2 (4 is also often used)
% 0017 %               INC=20;                             % set frame increment in samples
% 0018 %               NW=INC*OV;                          % DFT window length
% 0019 %               W=sqrt(hamming(NW,'periodic'));     % omit sqrt if OV=4
% 0020 %               W=W/sqrt(sum(W(1:INC:NW).^2));      % normalize window
% 0021 %               F=rfft(enframe(S,W,INC),NW,2);      % do STFT: one row per time frame, +ve frequencies only
% 0022 %               ... process frames ...
% 0023 %               X=overlapadd(irfft(F,NW,2),W,INC);  % reconstitute the time waveform (omit "X=" to plot waveform)
% 0024 %
% 0025 %  Inputs:   x    input signal
% 0026 %          win    window or window length in samples
% 0027 %          inc    frame increment in samples
% 0028 %            m    mode input:
% 0029 %                  'z'  zero pad to fill up final frame
% 0030 %                  'r'  reflect last few samples for final frame
% 0031 %                  'A'  calculate the t output as the centre of mass
% 0032 %                  'E'  calculate the t output as the centre of energy
% 0033 %
% 0034 % Outputs:   f    enframed data - one frame per row
% 0035 %            t    fractional time in samples at the centre of each frame
% 0036 %                 with the first sample being 1.
% 0037 %            w    window function used
% 0038 %

 nx=length(x(:));
 if nargin<2 || isempty(win)
     win=nx;
end
 if nargin<4 || isempty(m)
     m='';
 end
 nwin=length(win);
 if nwin == 1
     lw = win;
     w = ones(1,lw);
 else
     lw = nwin;
     w = win(:).';
 end
 if (nargin < 3) || isempty(inc)
     inc = lw;
 end
 nli=nx-lw+inc;
 nf = max(fix(nli/inc),0);   % number of full frames
 na=nli-inc*nf+(nf==0)*(lw-inc);       % number of samples left over
 fx=nargin>3 && (any(m=='z') || any(m=='r')) && na>0; % need an extra row
 f=zeros(nf+fx,lw);
 indf= inc*(0:(nf-1)).';
 inds = (1:lw);
 if fx
     f(1:nf,:) = x(indf(:,ones(1,lw))+inds(ones(nf,1),:));
     if any(m=='r')
         ix=1+mod(nf*inc:nf*inc+lw-1,2*nx);
         f(nf+1,:)=x(ix+(ix>nx).*(2*nx+1-2*ix));
     else
         f(nf+1,1:nx-nf*inc)=x(1+nf*inc:nx);
     end
     nf=size(f,1);
 else
     f(:) = x(indf(:,ones(1,lw))+inds(ones(nf,1),:));
 end
 if (nwin > 1)   % if we have a non-unity window
     f = f .* w(ones(nf,1),:);
 end
 if nargout>1
     if any(m=='E')
         t0=sum((1:lw).*w.^2)/sum(w.^2);
     elseif any(m=='A')
         t0=sum((1:lw).*w)/sum(w);
     else
         t0=(1+lw)/2;
     end
     t=t0+inc*(0:(nf-1)).';
 end