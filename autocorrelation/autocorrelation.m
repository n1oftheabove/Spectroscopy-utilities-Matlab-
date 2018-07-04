%%%%%%%% Autokorrelation mit BBO 1.0mm, theta 29.2°
%%% mit verschieden starken Filtern im Probe-Strahl
%% Datencharakterisierung
axisfontsize = 11;
%Messpunkte
N = 90000;
%Pulslaenge in fs
pulselength = 120;
axisfontsize = 15;

% % Zeitueberlapp bei stagepos (in mm)
% timezero = -18.17;
% % Daten gemessen bei stagepos (in mm)
% stagepos = -19.82;

% timeshift = 2* (abs(stagepos) - abs(timezero)) * 3.33;

% Einlesen
autocorr0 = load('autocorrelation 3500.dat','ascii');
autocorr10 = load('autocorrelation 11300 1,0.dat','ascii');

time0 = autocorr0(:,1);
data0 = autocorr0(:,2);
time10 = autocorr10(:,1);
data10 = autocorr10(:,2);


%% Daten nur von timemin bis timemax
% Für OD 0,0

[~,I1]=min(time0);
[~,I2]=max(time0);
if I2>I1
    a = 1;
else
    a=-1;
end

autocorrtimecut0ps = time0(I1:a:I2);
autocorrcut0 = data0(I1:a:I2);
%Bereich des Offsets definieren
zeroref_start0 = 1;
zeroref_end0 = 19290;


% Ref-Achse eichen (plot(autocorrcut0))
meanautocorr0 = mean(autocorrcut0(zeroref_start0:zeroref_end0));
autocorr0finx = -1*(autocorrcut0 - meanautocorr0);
autocorr0fin = autocorr0finx./max(autocorr0finx);

%% Daten nur von timemin bis timemax
% Für OD 1,0

[~,I1]=min(time10);
[~,I2]=max(time10);
if I2>I1
    a = 1;
else
    a=-1;
end

autocorrtimecut10 = time10(I1:a:I2);
autocorrcut10 = data10(I1:a:I2);

zeroref_start10 = 30930;
zeroref_end10 = 49250;

% Ref-Achse eichen (plot(autocorrcut0))
meanautocorr10 = mean(autocorrcut10(zeroref_start10:zeroref_end10));
autocorr10finx = -1*(autocorrcut10 - meanautocorr10);
autocorr10fin = autocorr10finx./max(autocorr10finx);

autocorrtimecut0 = autocorrtimecut0ps * 1000;


%% fd
a1=0.9807;
b1=1.968;
c1=100.1;
gaussy = a1*exp(-((autocorrtimecut0-b1)./c1).^2);
FWHM = 2*(c1*sqrt(-log(1./(2*a1)))+b1);
%%

figure(1);
s1=subplot(1,1,1);
diagramheight = 10; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
set(gca,'fontsize',axisfontsize,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
hold on
plot(autocorrtimecut0,autocorr0fin,'.r','markersize',0.1)
plot(autocorrtimecut0,gaussy,'.b','markersize',0.001)
%set(gca,'FontName','Times','FontSize',13);
hold off
axis xy;
axis(s1,[-700 700 -0.2 1.2]);
xlabel('Zeit in fs','Fontsize',axisfontsize);
ylabel('Normalisiertes Autokorrelationssignal','FontSize',axisfontsize);
grid
print(gcf,'autocorrelation','-depsc');


%%
figure(2);
plot(autocorrtimecut10,autocorr10fin,'.g','markersize',0.3)
xlabel('time in ps','fontsize',axisfontsize);
ylabel('Intensity (a.u.)','fontsize',axisfontsize);
title({'Autocorrelation signal','probe beam attenuated with ND 1.0'},'fontsize',axisfontsize-3);
grid
