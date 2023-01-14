clc;
clear all;
close all;
cd 5dB
[J P]=uigetfile('.wav','select the source wave file');
cd ..
[I Fs]=wavread(strcat(P,J));
figure,plot(I);grid on;xlabel('---Time');ylabel('---Amplitude');
title('Original speech signal');
%=====================================================
%===== Pre processing using MMSE ====================
[ss,gg,tt,ff,zo]=ssubmmsev(I,Fs);
figure,plot(ss);grid on;xlabel('---Time');ylabel('---Amplitude');
title('Pre-Processed speech signal');
%====================================================
%===================================================
[voiced, pitch_plot] = f_VOICED (ss,Fs,32e-3);
%=====================================================
Tw = 32; % analysis frame duration (ms) 
Ts = Tw/8; % analysis frame shift (ms)
[speechprocessed] = griffmet(ss,Fs,Tw,Ts,'G&L');
r=voiced.*speechprocessed(1:length(voiced));
T=r+speechprocessed(1:length(voiced));
figure,plot(T);grid on;xlabel('---Time');ylabel('---Amplitude');
title('Enhanced speech signal');
wavplay(I,Fs)
wavplay(T,Fs);
