%function_main of voiced/unvoiced detection
function [voiced, pitch_plot] = f_VOICED(x, fs, fsize);
f=1;
b=1;        %index no. of starting data point of current frame
frame_length = round(fs .* fsize);   %=number data points in each framesize of "x"
N= frame_length - 1;        %N+1 = frame length = number of data points in each framesize

%FRAME SEGMENTATION:
for b=1 : frame_length : (length(x) - frame_length),
    y1=x(b:b+N);     %"b+N" denotes the end point of current frame.
                %"y" denotes an array of the data points of the current 
                %frame
    y = filter([1 -.9378], 1, y1);  %pre-emphasis filter

    msf(b:(b + N)) = func_vd_msf (y);
    zc(b:(b + N)) = func_vd_zc (y);
    pitch_plot(b:(b + N)) = func_pitch (y,fs);
end

thresh_msf = (( (sum(msf)./length(msf)) - min(msf)) .* (0.67) ) + min(msf);
voiced_msf =  msf > thresh_msf;     %=1,0

thresh_zc = (( ( sum(zc)./length(zc) ) - min(zc) ) .*  (1.5) ) + min(zc);
voiced_zc = zc < thresh_zc;

thresh_pitch = (( (sum(pitch_plot)./length(pitch_plot)) - min(pitch_plot)) .* (0.5) ) + min(pitch_plot);
voiced_pitch =  pitch_plot > thresh_pitch;

for b=1:(length(x) - frame_length),
    if voiced_msf(b) .* voiced_pitch(b) .* voiced_zc(b) == 1,
%     if voiced_msf(b) + voiced_pitch(b) > 1,
        voiced(b) = 1;
    else
        voiced(b) = 0;
    end
end


