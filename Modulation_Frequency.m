
%% Comparative Acoustic Analysis of a Car Engine Condition Using Audio Signal Processing Techniques
%  Work developed by 
%% Muhammad Ali Saqib & Shahariar Ryehan
%  for the course
%% AUDIO PROCESSING AND SYNTHESIS  

clear;
clc; 
close all;

wavList=["normal_engine.wav","faulty_engine.wav"];  % (Please change the name of files to be analysed here)

cfg.winShort=1024;
cfg.winLong=4096;
cfg.overlap=0.90;
cfg.showCWT=true;

forw=wavList
assert(exist(w,'file')==2,['Filenotfound:',char(w)]);
end

featureTbl=table;
fork=1:numel(wavList)
feat=analyseOneFile(wavList(k),cfg);
featureTbl=[featureTbl;struct2table(feat)];
end

disp('');
disp('===================EXTRACTEDFEATURES===================');
disp(featureTbl);

functionfeat=analyseOneFile(wavFile,cfg)
[x,fs]=audioread(wavFile);
ifsize(x,2)>1,x=mean(x,2);end
x=detrend(x,'constant');
dur=numel(x)/fs;

fprintf('\n---%s(%.1fs@%dHz)---\n',char(wavFile),dur,fs);

mySTFT=@(sig,win,hop)stftOrSpect(sig,fs,win,hop);

winS=cfg.winShort;hopS=round((1-cfg.overlap)*winS);
winL=cfg.winLong;hopL=round((1-cfg.overlap)*winL);

[Ss,Fs,Ts]=mySTFT(x,winS,hopS);
[Sl,Fl,Tl]=mySTFT(x,winL,hopL);

env=abs(hilbert(x));
[Se,Fm,Tm]=mySTFT(env,winS,hopS);

figure('Name',['STFT/Envelope–',char(wavFile)],'NumberTitle','off');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
imagesc(Ts,Fs,mag2db(abs(Ss)));axisxy
title('STFT–Short');xlabel('Time(s)');ylabel('Hz');colorbar;caxis([-1100])

nexttile
imagesc(Tl,Fl,mag2db(abs(Sl)));axisxy
title('STFT–Long');xlabel('Time(s)');ylabel('Hz');colorbar;caxis([-1100])

nexttile
imagesc(Tm,Fm,mag2db(abs(Se)));axisxy
title('EnvelopeSTFT');xlabel('Time(s)');ylabel('ModHz');colorbar;caxis([-1100])

nexttile
ifexist('kurtosisMap','file')==2
kurtosisMap(env,fs);
title('Kurtosismap');
else
[Sz,Fz]=spectrogram(env,hamming(winS,'periodic'),winS-hopS,winS,fs,'yaxis');
SK=kurtosis(abs(Sz),0,2);
plot(SK,Fz,'LineWidth',1);set(gca,'YDir','reverse');gridon
xlabel('Spectralkurtosis');ylabel('Hz');
title('DIYspectralkurtosis');
end

ifcfg.showCWT&&exist('cwt','file')==2
figure('Name',['CWT–',char(wavFile)],'NumberTitle','off');
cwt(x,fs,'Wavelet','amor');
title(['CWTscalogram–',char(wavFile)]);
end

feat.File=string(wavFile);
feat.Duration_s=dur;
feat.RMS=rms(x);
feat.CrestFactor=max(abs(x))/feat.RMS;
feat.Skewness=skewness(x);
feat.Kurtosis=kurtosis(x);

[~,pkM]=max(abs(Se),[],1);
feat.ModFreq_Hz=median(Fm(pkM));

[~,pkC]=max(abs(Sl),[],1);
feat.CarrierFreq_Hz=median(Fl(pkC));
end

function[S,F,T]=stftOrSpect(x,fs,win,hop)
ifexist('stft','file')==2
[S,F,T]=stft(x,fs,'Window',hamming(win,'periodic'),'OverlapLength',win-hop,'FFTLength',win);
else
[S,F,T]=spectrogram(x,hamming(win,'periodic'),win-hop,win,fs,'yaxis');
end
end
