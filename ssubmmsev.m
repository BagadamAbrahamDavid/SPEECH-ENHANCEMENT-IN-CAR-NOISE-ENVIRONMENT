function [ss,gg,tt,ff,zo]=ssubmmsev(si,fsz,pp)
 %SSUBMMSE performs speech enhancement using mmse estimate of spectral amplitude or log amplitude [SS,ZO]=(S,FSZ,P)
%
 % Usage: y=ssubmmsev(x,fs);   % enhance the speech using default parameters
 %
% Inputs:
 %   si      input speech signal
 %   fsz     sample frequency in Hz
% 0009 %           Alternatively, the input state from a previous call (see below)
% 0010 %   pp      algorithm parameters [optional]
% 0011 %
% 0012 % Outputs:
% 0013 %   ss        output enhanced speech
% 0014 %   gg(t,f,i) selected time-frequency values (see pp.tf below)
% 0015 %   tt        centre of frames (in seconds)
% 0016 %   ff        centre of frequency bins (in Hz)
% 0017 %   zo        output state (or the 2nd argument if gg,tt,ff are omitted)
% 0018 %
% 0019 % The algorithm operation is controlled by a small number of parameters:
% 0020 %
% 0021 %        pp.of          % overlap factor = (fft length)/(frame increment) [2]
% 0022 %        pp.ti          % desired frame increment [0.016 seconds]
% 0023 %        pp.ri          % set to 1 to round ti to the nearest power of 2 samples [0]
% 0024 %        pp.ta          % time const for smoothing SNR estimate [0.396 seconds]
% 0025 %        pp.gx          % maximum posterior SNR as a power ratio [1000 = +30dB]
% 0026 %        pp.gn          % min posterior SNR as a power ratio when estimating prior SNR [1 = 0dB]
% 0027 %        pp.gz          % min posterior SNR as a power ratio [0.001 = -30dB]
% 0028 %        pp.xn          % minimum prior SNR [0]
% 0029 %        pp.xb          % bias compensation factor for prior SNR [1]
% 0030 %        pp.lg          % MMSE target: 0=amplitude, 1=log amplitude, 2=perceptual Bayes [1]
% 0031 %        pp.tn;         % smoothing time constant for noise estimation [0.5 s]
% 0032 %        pp.le;         % VAD threshold: log(p/(1-p)) where p is speech prob in a freq bin; use -Inf to prevent updating [0.15 (=>p=0.54)]
% 0033 %        pp.tx;         % initial noise interval [0.06 s]
% 0034 %        pp.ne          % noise estimation: 0=min statistics, 1=MMSE [0]
% 0035 %        pp.bt          % threshold for binary gain or -1 for continuous gain [-1]
% 0036 %        pp.mx          % input mixture gain [0]
% 0037 %        pp.rf          % round output signal to an exact number of frames [0]
% 0038 %        pp.tf          % selects time-frequency planes to output in the gg() variable ['g']
% 0039 %                           'i' = input power spectrum
% 0040 %                           'I' = input complex spectrum
% 0041 %                           'n' = noise power spectrum
% 0042 %                           'z' = "posterior" SNR (i.e. (S+N)/N )
% 0043 %                           'x' = "prior" SNR (i.e. S/N )
% 0044 %                           'g' = gain
% 0045 %                           'o' = output power spectrum
% 0046 %                           'O' = output complex spectrum
% 0047 %
% 0048 % The applied gain is mx+(1-mx)*optgain where optgain is calculated according to [1] or [2].
% 0049 % If pp.bt>=0 then optgain is first thresholded with pp.bt to produce a binary gain 0 or 1.
% 0050 %
% 0051 % The default parameters implement the original algorithm in [1,2].
% 0052 %
% 0053 % Several parameters relate to the estimation of xi, the so-called "prior SNR",
% 0054 %
% 0055 %             xi=max(a*pp.xb*xu+(1-a)*max(gami-1,pp.gn-1),pp.xn);
% 0056 %
% 0057 % This is estimated as a smoothed version of 1 less than gami, the "posterior SNR"
% 0058 % which is the noisy speech power divided by the noise power. This is
% 0059 % clipped to a min of (pp.gn-1), smoothed using a factor "a" which corresponds to a
% 0060 % time-constant of pp.ta and then clipped to a minimum of pp.xn. The
% 0061 % previous value is taken to be pp.xb*xu where xu is the ratio of the
% 0062 % estimated speech amplitude squared to the noise power.
% 0063 %
% 0064 % The noise estimation uses a VAD from equation (4) in [6] and recursively updates
% 0065 % the noise spectrum only in frames that are classified as noise-only.
% 0066 %
% 0067 % If convenient, you can call specsub in chunks of arbitrary size. Thus the following are equivalent:
% 0068 %
% 0069 %                   (a) y=ssubmmse(s,fs);
% 0070 %
% 0071 %                   (b) [y1,z]=ssubmmse(s(1:1000),fs);
% 0072 %                       [y2,z]=ssubmmse(s(1001:2000),z);
% 0073 %                       y3=ssubmmse(s(2001:end),z);
% 0074 %                       y=[y1; y2; y3];
% 0075 %
% 0076 % If the number of output arguments is either 2 or 5, the last partial frame of samples will
% 0077 % be retained for overlap adding with the output from the next call to ssubmmse().
% 0078 %
% 0079 % See also specsub() for an alternative gain function
% 0080 %
% 0081 % Refs:
% 0082 %    [1] Ephraim, Y. & Malah, D.
% 0083 %        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
% 0084 %        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
% 0085 %    [2] Ephraim, Y. & Malah, D.
% 0086 %        Speech enhancement using a minimum mean-square error log-spectral amplitude estimator
% 0087 %        IEEE Trans Acoustics Speech and Signal Processing, 33(2):443-445, Apr 1985
% 0088 %    [3] Rainer Martin.
% 0089 %        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
% 0090 %        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
% 0091 %    [4] O. Cappe.
% 0092 %        Elimination of the musical noise phenomenon with the ephraim and malah noise suppressor.
% 0093 %        IEEE Trans Speech Audio Processing, 2 (2): 345?349, Apr. 1994. doi: 10.1109/89.279283.
% 0094 %    [5] J. Erkelens, J. Jensen, and R. Heusdens.
% 0095 %        A data-driven approach to optimizing spectral speech enhancement methods for various error criteria.
% 0096 %        Speech Communication, 49: 530?541, 2007. doi: 10.1016/j.specom.2006.06.012.
% 0097 %    [6] J. Sohn, N. S. Kim, and W. Sung.
% 0098 %        A statistical model-based voice activity detection.
% 0099 %        IEEE Signal Processing Lett., 6 (1): 1?3, 1999. doi: 10.1109/97.736233.
% 0100 %    [7] Loizou, P.
% 0101 %        Speech enhancement based on perceptually motivated Bayesian estimators of the speech magnitude spectrum.
% 0102 %        IEEE Trans. Speech and Audio Processing, 13(5), 857-869, 2005.
% 0103 

 persistent kk cc 
 if ~numel(kk)
     kk=sqrt(2*pi);      % sqrt(8)*Gamma(1.5) - required constant
     cc=sqrt(2/pi);      %sqrt(2)/Gamma(0.5)
 end
 if numel(si)>length(si)
     error('Input speech signal must be a vector not a matrix');
 end
 if isstruct(fsz)
     fs=fsz.fs;
     qq=fsz.qq;
     qp=fsz.qp;
     ze=fsz.ze;
     s=zeros(length(fsz.si)+length(si(:)),1); % allocate space for speech
     s(1:length(fsz.si))=fsz.si;
     s(length(fsz.si)+1:end)=si(:);
 else
     fs=fsz;     % sample frequency
     s=si(:);
     % default algorithm constants
     
     qq.of=2;        % overlap factor = (fft length)/(frame increment)
     qq.ti=16e-3;    % desired frame increment (16 ms)
     qq.ri=0;        % round ni to the nearest power of 2
     qq.ta=0.396;    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [1]
     qq.gx=1000;     % maximum posterior SNR = 30dB
     qq.gn=1;        % min posterior SNR as a power ratio when estimating prior SNR [1]
     qq.gz=0.001;    % min posterior SNR as a power ratio [0.001 = -30dB]
     qq.xn=0;        % minimum prior SNR = -Inf dB
     qq.xb=1;        % bias compensation factor for prior SNR [1]
     qq.lg=1;        % use log-domain estimator by default
     qq.ne=0;        % noise estimation: 0=min statistics, 1=MMSE [0]
     qq.bt=-1;       % suppress binary masking
     qq.mx=0;        % no input mixing
     qq.tf='g';      % output the gain time-frequency plane by default
     qq.rf=0;
     qq.tn=0.5;     % smoothing constant for noise estimation [500 ms]
     qq.le=0.15;    % VAD threshold; use -Inf to prevent updating
     qq.tx=0.06;    % initial noise interval [60 ms]
     if nargin>=3 && ~isempty(pp)
         qp=pp;      % save for estnoisem call
         qqn=fieldnames(qq);
         for i=1:length(qqn)
         if isfield(pp,qqn{i})
                 qq.(qqn{i})=pp.(qqn{i});
             end
         end
     else
         qp=struct;  % make an empty structure
     end
 end
 % derived algorithm constants
 if qq.ri
   ni=pow2(nextpow2(ti*fs*sqrt(0.5)));
 else
     ni=round(qq.ti*fs);    % frame increment in samples
 end
 tinc=ni/fs;         % true frame increment time
 a=exp(-tinc/qq.ta); % SNR smoothing coefficient
 gx=qq.gx;           % max posterior SNR as a power ratio
 gz=qq.gz;           % min posterior SNR as a power ratio
 xn=qq.xn;           % floor for prior SNR, xi
 ne=qq.ne;           % noise estimation: 0=min statistics, 1=MMSE [0]
 gn1=max(qq.gn-1,0); % floor for posterior SNR when estimating prior SNR
 le=qq.le;
 xb=qq.xb;
 tf=qq.tf;
 rf=qq.rf || nargout==2 || nargout==5;            % round down to an exact number of frames
 nd=max(1,round(qq.tx/tinc)); % number of frames for initial noise estimate
 an=exp(-tinc/qq.tn); % Noise spectrum smoothing coefficient
 
%% calculate power spectrum in frames
 
 no=round(qq.of);                    % integer overlap factor
 nf=ni*no;                           % fft length
 w=sqrt(hamming(nf+1))'; w(end)=[];  % for now always use sqrt hamming window
 w=w/sqrt(sum(w(1:ni:nf).^2));       % normalize to give overall gain of 1
 if rf>0
    rfm='';                         % truncated input to an exact number of frames
 else
     rfm='r';
 end
 [y,tt]=enframe(s,w,ni,rfm);
 tt=tt/fs;                           % frame times
 yf=rfft(y,nf,2);
 yp=yf.*conj(yf);                    % power spectrum of input speech
 [nr,nf2]=size(yp);                  % number of frames
 ff=(0:nf2-1)*fs/nf;
 if isstruct(fsz)
    ndp=fsz.ndp;
     dpi=fsz.dpi;
     ssv=fsz.ssv;
     xu=fsz.xu;                      % saved unsmoothed SNR
 else
     dpi=zeros(1,nf2);   % noise estimate
     ndp=0;              % noise estimate based on ndp frames
     ssv=zeros(ni*(no-1),1);            % dummy saved overlap
     xu=1;                           % dummy unsmoothed SNR from previous frame
 end
 if ~nr                                 % no data frames
     ss=[];
     gg=[];
 else
     if ndp<nd
         ndx=min(nr,nd-ndp);         % number of frames to use
         dpi=ndp/(ndp+ndx)*dpi+sum(yp(1:ndx,:),1)/(ndp+ndx);
         ndp=ndp+ndx;
     end
     g=zeros(nr,nf2);                % create space for gain matrix
     x=zeros(nr,nf2);                % create space for prior SNR
     dp=zeros(nr,nf2);               % create space for noise power spectrum estimate
     switch qq.lg
         case 0                      % use amplitude domain estimator from [1]
             for i=1:nr
                 ypi=yp(i,:);
                 gami=max(min(ypi./dpi,gx),gz);     % gamma = posterior SNR
                 xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
                 if sum(gami.*xi./(1+xi)-log(1+xi))<le*nf2 % noise frame
                     dpi=dpi*an+(1-an)*ypi;
                 end
                 dp(i,:)=dpi;  % only required if noise estimate output is requested
                 v=0.5*xi.*gami./(1+xi);    % note that this is 0.5*vk in [1]
                 gi=(0.277+2*v)./gami;     % accurate to 0.02 dB for v>0.5
                 mv=v<0.5;
                 if any(mv)
                     vmv=v(mv);
                     gi(mv)=kk*sqrt(vmv).*((0.5+vmv).*besseli(0,vmv)+vmv.*besseli(1,vmv))./(gami(mv).*exp(vmv));
                 end
                 g(i,:)=gi;              % save gain for later
                 x(i,:)=xi;              % save prior SNR
                 xu=gami.*gi.^2;         % unsmoothed prior SNR
             end
         case 2                          % perceptually motivated estimator from [7]
             for i=1:nr
                 ypi=yp(i,:);
                 gami=max(min(ypi./dpi,gx),gz);     % gamma = posterior SNR
                xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
                 if sum(gami.*xi./(1+xi)-log(1+xi))<le*nf2 % noise frame
                 dpi=dpi*an+(1-an)*ypi;
            end
              v=0.5*xi.*gami./(1+xi);    % note that this is 0.5*vk in [7]
              gi=cc*sqrt(v).*exp(v)./(gami.*besseli(0,v));
                g(i,:)=gi;              % save gain for later
                x(i,:)=xi;              % save prior SNR
                xu=gami.*gi.^2;         % unsmoothed prior SNR
             end
      otherwise                       % use log domain estimator from [2]
             for i=1:nr
                 ypi=yp(i,:);
                 gami=max(min(ypi./dpi,gx),gz);     % gamma = posterior SNR
                 xi=max(a*xb*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
                 xir=xi./(1+xi);
                 if sum(gami.*xir-log(1+xi))<le*nf2 % noise frame
                     dpi=dpi*an+(1-an)*ypi;
                 end
                 gi=xir.*exp(0.5*expint(xir.*gami));
                g(i,:)=gi;                 % save gain for later
                 x(i,:)=xi;              % save prior SNR
                xu=gami.*gi.^2;         % unsmoothed prior SNR
          end
    end
    if qq.bt>=0
      g=g>qq.bt;
  end
     g=qq.mx+(1-qq.mx)*g;                    % mix in some of the input
     se=(irfft((yf.*g).',nf).').*repmat(w,nr,1);     % inverse dft and apply output window
     ss=zeros(ni*(nr+no-1),no);                      % space for overlapped output speech
     ss(1:ni*(no-1),end)=ssv;
     for i=1:no
        nm=nf*(1+floor((nr-i)/no));         % number of samples in this set
         ss(1+(i-1)*ni:nm+(i-1)*ni,i)=reshape(se(i:no:nr,:)',nm,1);
     end
    ss=sum(ss,2);
    if nargout>2 && ~isempty(tf)
         gg=zeros(nr,nf2,length(tf));  % make space
         for i=1:length(tf)
             switch tf(i)
                 case 'i'            % 'i' = input power spectrum
                     gg(:,:,i)=yp;
                case 'I'            % 'i' = input power spectrum
                     gg(:,:,i)=yf;
                 case 'n'            % 'n' = noise power spectrum
                     gg(:,:,i)=dp;
                 case 'z'            % 'z' = posterior SNR (i.e. (S+N)/N )
                     gg(:,:,i)=gam;
                 case 'x'            % 'x' = prior SNR
                     gg(:,:,i)=x;
                 case 'g'            % 'g' = gain
                     gg(:,:,i)=g;
                 case 'o'            % 'o' = output power spectrum
                     gg(:,:,i)=yp.*g.^2;
                 case 'O'            % 'o' = output power spectrum
                     gg(:,:,i)=yf.*g;
             end
         end
     end
 end
 if nargout==2 || nargout==5
     if nr
         zo.ssv=ss(end-ni*(no-1)+1:end);    % save the output tail for next time
         ss(end-ni*(no-1)+1:end)=[];        % only output the frames that are completed
     else
         zo.ssv=ssv;  %
     end
     zo.si=s(length(ss)+1:end);      % save the tail end of the input speech signal
     zo.fs=fs;                       % save sample frequency
     zo.qq=qq;                       % save local parameters
     zo.qp=qp;                       % save estnoisem parameters
     zo.xu=xu;
     zo.dpi=dpi;
     zo.ndp=ndp;
     if nargout==2
         gg=zo;                      % 2nd of two arguments is zo
     end
 elseif rf==0
     ss=ss(1:length(s));             % trim to the correct length if not an exact number of frames
 end
 if ~nargout && nr>0
     ffax=ff/1000;
     ax=zeros(4,1);
     ax(1)=subplot(223);
     imagesc(tt,ffax,20*log10(g)');
     colorbar;
     axis('xy');
     title(sprintf('Filter Gain (dB): ta=%.2g',qq.ta));
     xlabel('Time (s)');
     ylabel('Frequency (kHz)');
     
     ax(2)=subplot(222);
     imagesc(tt,ffax,10*log10(yp)');
     colorbar;
     axis('xy');
     title('Noisy Speech (dB)');
     xlabel('Time (s)');
     ylabel('Frequency (kHz)');
     
     ax(3)=subplot(224);
     imagesc(tt,ffax,10*log10(yp.*g.^2)');
     colorbar;
     axis('xy');
     title('Enhanced Speech (dB)');
     xlabel('Time (s)');
     ylabel('Frequency (kHz)');
     
     ax(4)=subplot(221);
     imagesc(tt,ffax,10*log10(dp)');
     colorbar;
     axis('xy');
    title('Noise Estimate (dB)');
     xlabel('Time (s)');
     ylabel('Frequency (kHz)');
     linkaxes(ax);
end