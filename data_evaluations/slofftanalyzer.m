%this function demands a set of (non double) reflectivity and aequidistant
%delayvalues ('reflectivity' and 'delay'), a lower index 'k' and an upper
%index 'l', cuts the data in that index range, produces a smoothed one with
%a given 'smoothrange', subtracts the unsmoothed with the smoothed one,
%applies a window function and outputs the FFT of this ('fftcut' over
%'fftcut'). Also for the FFT provide information for the pulsewidth
%('pulsew') remark: Window-function 'window_1.m' necessary
function [ fftcutfreq, fftcut, cutdelay, rdy2windowref ] = slofftanalyzer(delay, reflectivity, k, l, smoothrange, pulsew)
cutdelay = delay(k:l);
cutreflectivity =  reflectivity(k:l);

smoothedreflectivity = smooth(cutreflectivity,smoothrange);
rdy2windowref = cutreflectivity - smoothedreflectivity;

win_width = 0.99*abs(max(cutdelay)-min(cutdelay));
win_cent = mean(cutdelay);
win_edge = 0.02;
window = window_1(win_width, win_edge, win_cent,cutdelay);
windwd_corcut = window.*rdy2windowref;

if max(delay).*min(delay) < 0
    delayrange = abs(max(delay))+abs(min(delay));
else
    delayrange = abs(max(delay))-abs(min(delay));
end
samples= length(delay);

T = delayrange./samples;      %time intervall between two measured points
Fs = 1/T;       % Sampling Frequency in THz
      
NFFT = length(cutdelay);
Y = fft(windwd_corcut, NFFT)/NFFT;
f = Fs/2 * linspace(0,1,floor(NFFT/2)+1);
Ampdens = 2*abs(Y(1:floor(NFFT/2)+1))';
%erster Index in frequenzvektor, ab dem Frequenz gr��er als Pulsbreite
%klein ist
fmax = find(f >(1000*1./pulsew),1,'first');
fftcutfreq= f(1:fmax);
fftcut= Ampdens(1:fmax);


% Fig1=figure;
% figure(Fig1);
% diagramheight = 10; %in cm in gedrucktem Dokument
% diagramwidth = 14; %in cm in gedrucktem Dokumen
% set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
% set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
% set(gcf,'PaperSize',[diagramwidth diagramheight]);
% set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
% plot(cutdelay,cutreflectivity,'.b','markersize',1);
% %title({'Sample STO/SRO on STO, 30 DL'},'fontsize',axisfontsize);
% xlabel('Delay in ps','fontsize',15);
% ylabel('\DeltaR/R_0','fontsize',15);
% grid
% 
% Fig2=figure;
% figure(Fig2);
% diagramheight = 10; %in cm in gedrucktem Dokument
% diagramwidth = 14; %in cm in gedrucktem Dokumen
% set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
% set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
% set(gcf,'PaperSize',[diagramwidth diagramheight]);
% set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
% plot(cutdelay,rdy2windowref,'.b','markersize',1);
% %title({'Sample STO/SRO on STO, 30 DL'},'fontsize',axisfontsize);
% xlabel('Delay in ps','fontsize',15);
% ylabel('\DeltaR/R_0','fontsize',15);
% grid
% 
% Fig3=figure;
% figure(Fig3);
% diagramheight = 10; %in cm in gedrucktem Dokument
% diagramwidth = 14; %in cm in gedrucktem Dokumen
% set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
% set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
% set(gcf,'PaperSize',[diagramwidth diagramheight]);
% set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
% hold on
% plot(cutdelay, windwd_corcut,'.b','markersize',1);
% plot(cutdelay, window*2*1e-6,'.r','markersize',1);
% hold off
% %title({'Sample STO/SRO on STO, 30 DL'},'fontsize',axisfontsize);
% xlabel('Delay in ps','fontsize',15);
% ylabel('\DeltaR/R_0','fontsize',15);
% grid
% 
% 
% Fig4=figure;
% figure(Fig4);
% diagramheight = 10; %in cm in gedrucktem Dokument
% diagramwidth = 14; %in cm in gedrucktem Dokumen
% set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
% set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
% set(gcf,'PaperSize',[diagramwidth diagramheight]);
% set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
% plot(fftcutfreq(1:(find(fftcutfreq >1.9,1,'first'))), fftcut(1:(find(fftcutfreq >1.9,1,'first'))),'r','markersize',5,'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k');
% %title({'Sample STO/SRO on STO, 30 DL FFT'},'fontsize',axisfontsize);
% xlabel('Frequenz in THz','fontsize',15);
% ylabel('Amplitudendichte','fontsize',15);
% grid

% figure(1);
% subplot(3,1,1)
% plot(cutdelay,cutreflectivity,'.b','markersize',1);
% title('cut area of interest (AOE)','FontWeight','bold');
% subplot(3,1,2)
% plot(cutdelay,rdy2windowref,'.b','markersize',1);
% title('cut AOE, hard smoothed substracted','FontWeight','bold');
% subplot(3,1,3)
% hold on
% plot(cutdelay, windwd_corcut,'.b','markersize',1);
% plot(cutdelay, window*2*1e-6,'.r','markersize',1);
% title('cut AOE, hard smoothed substr, multiplied w. Window','FontWeight','bold');
% hold off
% figure(2);
% plot(fftcutfreq,fftcut,'r','markersize',1,'Marker','.');
% title('FFT Amplitude Spectrum','FontWeight','bold');
% xlabel('Frequency (THz)')
% ylabel('|Y(f)|')
% grid

