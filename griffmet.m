% G&L method
function [y] = griffmet(s, fs, Tw, Ts, stype)

    %______________________________________________________________________________________________________________________________
    if(nargin<4), error(sprintf('Not enough input arguments. Type "help %s" for usage help.', mfilename)); end;
    if(nargin==4), stype='GRIFFIN & LIM','G&L'; end;

    %______________________________________________________________________________________________________________________________
    s = s(:).'-mean(s);                              % make sure input signal is in row form and zero-mean
    Nw = round(fs*Tw*0.001);                         % frame duration (in samples)
    Ns = round(fs*Ts*0.001);                         % frame shift (in samples)
    nfft = 2^nextpow2(2*Nw);                         % FFT analysis length

    winfunc = @(L,S)(sqrt(2*S/(L*(2*0.5^2+(-0.5)^2)))*(0.5-0.5*cos((2*pi*((0:L-1)+0.5))/L)));
    w = winfunc(Nw, Ns);                             % Griffin & Lim's modified Hanning window 
    %w = hamming(Nw,'periodic').';                    % periodic Hamming window

    D = mod(length(s), Ns);                          % we will add Nw-D zeros to the end
    G = (ceil(Nw/Ns)-1)*Ns;                          % we will add G zeros to the beginning
    s = [zeros(1,G) s zeros(1,Nw-D)];                % zero pad signal to allow an integer number of segments
    L = length(s);                                   % length of the signal for processing (after padding)
    M = ((L-Nw)/Ns)+1;                               % number of overlapped segments

    %______________________________________________________________________________________________________________________________
    indf = Ns*[0:(M-1)].';                           % frame indices
    inds = [1:Nw];                                   % sample indices in each frame
    refs = indf(:,ones(1,Nw)) + inds(ones(M,1),:);   % absolute sample indices for each frame
    frames = s(refs);                                % split into overlapped frames using indexing (frames as rows)
    frames_w = frames.*repmat(w, M, 1);              % apply analysis window
    S = fft(frames_w, nfft, 2);                      % perform short-time Fourier transform (STFT) analysis
    MAG = abs(S);                                    % compute STFT magnitude spectra (across each row/frame)
    PHA = angle(S);                                  % compute STFT phase spectra (across each row/frame)
    MSTFS = MAG.*exp(j*PHA);                         % recombine magnitude and phase to produce modified STFT
    %______________________________________________________________________________________________________________________________
    x = real(ifft(MSTFS, nfft, 2));                 % perform inverse STFT analysis
    x = x(:, 1:Nw);                                 % discard FFT padding from frames

    switch(upper(stype)) 
    case {'GRIFFIN & LIM','G&L'}                     % Griffin & Lim's method

        x = x .* w(ones(M,1),:);                     % apply synthesis window (Griffin & Lim's method)
        y = zeros(1, L); 
        for i = 1:M, y(refs(i,:)) = y(refs(i,:)) + x(i,:); end; % overlap-add processed frames
        wsum2 = zeros(1, L);
        for i = 1:M, wsum2(refs(i,:)) = wsum2(refs(i,:)) + w.^2; end; % overlap-add squared window samples
        y = y./wsum2;                                % divide out squared and summed-up analysis windows

    otherwise, error(sprintf('%s: synthesis type not supported.', stype));
    end

    y = y(G+1:L-(Nw-D));                         % remove the padding


