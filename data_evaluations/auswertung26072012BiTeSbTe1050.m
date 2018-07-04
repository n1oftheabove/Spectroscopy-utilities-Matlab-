%% %% BiTe (10Angstroem) / SbTe (50Angstroem) on GaAs
%% Input Values
%Fuer welche Wellenlaengen liegen Messungen vor?
wvl = [795];
wvllist = {'795 nm'};
%Records?
records = 10000 ;
%Messwerte pro Umlauf
timesteps = 90000;
%Zeitueberlapp bei stagepos (in mm)
timezero = -18.140;
%Daten gemessen bei stagepos (in mm)
stagepos = -18.64;
timeshift = 2* (abs(stagepos) - abs(timezero)) * (1e9./physconst('lightspeed'));

%Indexshift
timeindexshift = -21;

wvlsize = length(wvl);
pulselength = 0.120;

fftmarkersize = 5;

%% Erstellen einer Matrix mit den existierenden Datei-Namen
 for i = 1:wvlsize,
 databack_name{i} = sprintf('%dnm %d back.dat', wvl(i), records);
 end;
 
for i = 1:wvlsize,
data_name{i} = sprintf('%dnm %d.dat', wvl(i), records);
end;

back = zeros(timesteps,wvlsize);
data = zeros(timesteps,wvlsize);
temp = zeros(timesteps,2);
time = zeros(timesteps,1);

%% Einlesen und Variablendefinieren der Messmatrizen fuer background
% echte Messung

 for i=1:wvlsize
     temp = load(databack_name{i},'-ascii');
     back(:,i) = temp(:,2);
     time = temp(:,1);
 end

temp = zeros(timesteps,2);

for i=1:wvlsize
    temp = load(data_name{i},'-ascii');
    data(:,i) = temp(:,2);
end

% Korrigieren der Messwerte um die background-Messungen
cor = data - back;

%'Umklappen' der Messwerte - spaetere Zeiten kommen an den Anfang

% timeshifted = circshift(time,[timeindexshift 0]);
% time2 = timeshifted(1:end-abs(timeindexshift),:);
% cor2 = cor(1:end-abs(timeindexshift),:); 


%% Daten nur von timemin bis timemax
[timecut,cor2_1] = timecutter(time,cor(:,1));
cutrefmatrix = zeros(length(cor2_1),wvlsize);
for i = 1:wvlsize
    [~,cutrefmatrix(:,i)]=timecutter(time,cor(:,i));
end
%doppelte Zeitwerte rauskicken
[time_1,ref_1 ] = aequidisttime(timecut,cutrefmatrix(:,1));
refmatrix = zeros(length(ref_1),wvlsize);
for i = 1:wvlsize
    [~,refmatrix(:,i)]=aequidisttime(timecut,cutrefmatrix(:,i));
end
%ref Achse eichen 1:11260
zeroedrefmatrix = zeros(length(ref_1),wvlsize);
for i = 1:wvlsize
    zeroedrefmatrix(:,i) = refmatrix(:,i)-mean(refmatrix(1:11260,i));
end
% t-Achse eichen
time_fixed = time_1 + timeshift;
% plotten:
cc = flipud(autumn(wvlsize));
figure(2)
hold on
h = zeros(3,wvlsize);
for i = 1:wvlsize;
h(:,i)=plot(time_fixed,zeroedrefmatrix(:,i),'color',cc(i,:));
set(h(:,i), 'Color', cc(i,:))
end
hold off
legend(h(1,:),wvllist)
grid
% smoothen mit Pulslänge
delta = abs(max(time_fixed)-min(time_fixed))./length(time_fixed);
smoothedrefmatrix = zeros(length(zeroedrefmatrix(:,1)),wvlsize);
for i = 1:wvlsize;
    smoothedrefmatrix(:,i) = smooth(zeroedrefmatrix(:,i), pulselength./delta);
end

%% Interessanten Bereich ausschneiden und FFT davon machen (plot(zeroedrefmatrix))
% auch von mit Pulslänge geglättetem Spektrum
slofftanalyzerstart = 14200;
slofftanalyzerend = 48170;
fftsmspan = 4000;
[fftxx, fftyy, sloxx, sloyy] = slofftanalyzer(time_fixed, smoothedrefmatrix(:,1), slofftanalyzerstart, slofftanalyzerend, fftsmspan, 120);
fftxxmatrix = zeros(length(fftxx),wvlsize);
fftyymatrix = zeros(length(fftyy),wvlsize);
sloxxmatrix = zeros(length(sloxx),wvlsize);
sloyymatrix = zeros(length(sloyy),wvlsize);
for i=1:wvlsize
[fftxxmatrix(:,i), fftyymatrix(:,i), sloxxmatrix(:,i), sloyymatrix(:,i)] = slofftanalyzer(time_fixed, smoothedrefmatrix(:,i), slofftanalyzerstart, slofftanalyzerend, fftsmspan, 120);
end

% %% Alle FFTs plotte
% cc = flipud(autumn(wvlsize));
% figure(2)
% hold on
% h = zeros(3,wvlsize);
% for i = 1:wvlsize;
% h(:,i)=plot(fftxxmatrix(:,i),fftyymatrix(:,i),'color',cc(i,:));
% set(h(:,i), 'Color', cc(i,:))
% end
% hold off
% legend(h(1,:),wvllist)
% grid

%% Transiente plotten

smoothedrefmatrixBiTeSbTe1050 = smoothedrefmatrix(:,i);
time_fixedBiTeSbTe1050 = time_fixed;

cc = flipud(autumn(wvlsize));
fig=figure(1);
diagramheight = 10; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
set(gca,'fontsize',8,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
hold on
h = zeros(3,wvlsize);
for i = 1:wvlsize;
h(:,i)=plot(time_fixed,smoothedrefmatrix(:,i),'color',cc(i,:),'Markersize',3);
set(h(:,i), 'Color', cc(i,:))
end
hold off
%legend(h(1,:),wvllist)
grid
%title({'[BiTe (10angstrom)/SbTe (50angstrom)] on GaAs','10000 records'})
xlabel('Delay \tau in ps','FontSize',11)
ylabel('\DeltaR/R_0','FontSize',11)
print(fig, '26072012BiTeSbTe(10_50)10000', '-depsc');


%% Nur SLO plotte
cc = flipud(autumn(wvlsize));
fig=figure(2);
diagramheight = 10; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
hold on
h = zeros(3,wvlsize);
for i = 1:wvlsize;
h(:,i)=plot(sloxxmatrix(:,i),sloyymatrix(:,i),'color',cc(i,:));
set(h(:,i), 'Color', cc(i,:))
end
hold off
%legend(h(1,:),wvllist)
grid
%title({'[BiTe (10angstrom)/SbTe (50angstrom)] on GaAs','Oscillations only 10000 records'})
xlabel('Zeit in ps')
ylabel('\DeltaR/R_0')
print(fig, '26072012BiTeSbTe(10_50)10000slo', '-depsc');

%% Beide in eines neu

cutrefsmnan = nan(length(smoothedrefmatrix(:,1)),1);
cutdelaysmnan = nan(length(time_fixed),1);
cutrefsmnan(slofftanalyzerstart:slofftanalyzerend,1) = sloyymatrix(:,1);
cutdelaysmnan(slofftanalyzerstart:slofftanalyzerend,1) = sloxxmatrix(:,1);

h=figure;
diagramheight = 14; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
Ticksh1 = [-8e-4 -6e-4 -4e-4 -2e-4 0];
Ticksh2 = [-6e-5 -4e-5 -2e-5 0 2e-5 4e-5];
%Xticks=
Achsenh1 = [-2 12 -8e-4 1e-4];
Achsenh2 = [-2 12 -6e-5 6e-5];
Labelsh2 = [-0.6 -0.4 -0.2 0 0.2 0.4];
set(h,'units','centimeters');
%set(h,'position',[5 5 12 12]);
set(h,'PaperSize',[14 14]);
set(h,'PaperPositionMode','auto');

h2=subplot(2,2,1);
plot(time_fixed, cutrefsmnan,'.c','markersize',.5)
axis xy;
box on;
axis(h2,Achsenh2);
set(gca,'YTick',Ticksh2);
%set(gca,'XTick',Xticks);
set(gca,'YTickLabel',Labelsh2);
set(gca,'FontSize',8,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
grid

h1=subplot(2,2,2);
plot(time_fixed, smoothedrefmatrix(:,1),'.c','markersize',0.5);
set(gca,'YTick',Ticksh1);
set(gca,'XTickLabel',[]);
%set(gca,'XTick',Xticks);
set(gca,'FontSize',8,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
axis(h1,Achsenh1);
grid

h3=subplot(2,2,[3 4]);
box on
cc = flipud(autumn(wvlsize));
h = zeros(3,wvlsize);
hold on
for i = 1:wvlsize;
h(:,i)=plot(fftxxmatrix(:,1),fftyymatrix(:,i),'color',cc(i,:),'markersize',5,'Marker','.','MarkerEdgecolor','k');
set(h(:,i),'Color',cc(i,:));
end
hold off
set(gca,'FontSize',8,'XTick',[0 1 2 3 4 5]);
xlabel('\nu in THz','FontSize',11)
ylabel('|Y(\nu)|','FontSize',11)
axis xy;
axis(h3,[0 6 0 1e-5]);

set(h2,'position',[0.15 0.1 0.75 0.5]);
set(h1,'position',[0.15 0.6 0.75 0.24]);
set(h3,'position',[0.55 0.52 0.3 0.15]);
xlabel(h2,'Delay \tau in ps','FontSize',11);
ylabel(h1,'\DeltaR/R_0','FontSize',11);
ylabel(h2,'\DeltaR/R_0','FontSize',11);
grid
text(2.5, 0.7e-5,'\bfFFT','FontSize',13, 'BackgroundColor',[1 1 1]);
U = uipanel('backgroundcolor','w','pos',get(gca,'outerpos'),'BorderType','etchedout','ForeGroundColor',[1 1 1]);
CH = get(gcf,'children');
set(gcf,'children',CH([2 1 3 4]));
print(gcf, '26072012BiTeSbTe(10_50)10000doubleplot', '-depsc');

%% 
BiTeSbTeSLOxx = sloxxmatrix(:,1);
BiTeSbTeSLOyy = sloyymatrix(:,1);
BiTeSbTeSLOyy(1:514) = NaN;
BiTeSbTeSLOxxx = sloxxmatrix(514:end,1);
BiTeSbTeSLOyyy = sloyymatrix(514:end,1);

%% fit selber machen
A=-7.198e-5;
Gamma=0.5481;
Omega=12.67;
beta=0.0188;
phi=0.04069;
fig7=figure(7);
s1=subplot(1,1,1);
diagramheight = 7; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
set(gca,'fontsize',8,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');

y = A.*exp(-Gamma.*BiTeSbTeSLOxxx).*cos((Omega+beta.*BiTeSbTeSLOxxx).*BiTeSbTeSLOxxx+phi);

hold on
plot(BiTeSbTeSLOxxx, BiTeSbTeSLOyyy,'-r','MarkerSize',1)
plot(BiTeSbTeSLOxxx, y,'-b','MarkerSize',1)
%title({'[Bi_2Te_3 (10Angstroem) / Sb_2Te_3 (50Angstroem)] x 175 on GaAs (111)','\Phi_{pump} \approx 110 \muJ / cm^2'},'FontSize',12)
axis xy;
axis(s1,[0 12 -5e-5 6e-5])
xlabel('Delay \tau in ps','FontSize',11)
ylabel('\DeltaR / R_0','FontSize',11)
hleg1=legend('Gemessen','Fit-Funktion, s. Text');
legend('boxoff');
hLine = findobj(hleg1,'type','line');
set(hLine,'LineWidth',1.5);
hold off
grid
print(fig7, '26072012BiTeSbTe(10_50)10000plotwithfit', '-depsc');
%% FFT

fig=figure(6);
diagramheight = 10; %in cm in gedrucktem Dokument
diagramwidth = 14; %in cm in gedrucktem Dokumen
set(gcf,'Units','centimeter','Position',[0 0 diagramwidth diagramheight]);
set(gcf, 'PaperPositionMode','auto','PaperUnits','centimeter');
set(gcf,'PaperSize',[diagramwidth diagramheight]);
set(gca,'fontsize',13,'fontname','Helvetica', 'XMinorTick','on','YMinorTick','on');
s1=subplot(1,1,1);
cc = flipud(autumn(wvlsize));
h = zeros(3,wvlsize);
hold on
for i = 1:wvlsize;
h(:,i)=plot(fftxxmatrix(:,1),fftyymatrix(:,i),'color',cc(i,:),'markersize',fftmarkersize,'Marker','.','MarkerEdgecolor','k');
set(h(:,i),'Color',cc(i,:));
end
hold off

xlabel('\nu in THz')
ylabel('Amplitudendichte |Y|')
axis xy;
axis(s1,[0 8 0 1e-5]);
grid
print(fig, '26072012BiTeSbTe(10_50)10000fft', '-depsc');
%% Alles in einen Plot
fig=figure(10);
%annotation('textbox', [0 0.9 1 0.1], ...
%    'String',{'[Bi_2Te_3 (10Angstroem) / Sb_2Te_3 (50Angstroem)] x 175 + 2000Angstroem Bi_2Te_3 on GaAs (111)' ,'\Phi_{pump} \approx 123 \muJ / cm^2'},'FontSize',13,'EdgeColor', 'none','HorizontalAlignment', 'center')
%title({'[Bi_2Te_3 (20Angstroem) / Sb_2Te_3 (30Angstroem)] x 170 +','2000Angstroem Bi_2Te_3 on GaAs (111)','\Phi_{pump} \approx 110 \muJ / cm^2'},'FontSize',12)
set(gcf,'Units','pixel','Position',[0 0 1024 700]);
subplot(2,2,1)
cc = flipud(autumn(wvlsize));
hold on
h = zeros(3,wvlsize);
for i = 1:wvlsize;
h(:,i)=plot(time_fixed,smoothedrefmatrix(:,i),'color',cc(i,:),'Markersize',3);
set(h(:,i), 'Color', cc(i,:))
end
hold off
%legend(h(1,:),wvllist)
grid
xlabel('Zeit in ps','FontSize',13)
ylabel('\DeltaR/R_0','FontSize',13)

subplot(2,2,2)
cc = flipud(autumn(wvlsize));
hold on
h = zeros(3,wvlsize);
for i = 1:wvlsize;
h(:,i)=plot(sloxxmatrix(:,i),sloyymatrix(:,i),'color',cc(i,:));
set(h(:,i), 'Color', cc(i,:))
end
hold off
grid
xlabel('Zeit in ps','FontSize',13)
ylabel('\DeltaR/R_0','FontSize',13)

subplot(2,2,3)
hold on
plot(BiTeSbTeSLOxxx, BiTeSbTeSLOyyy,'.b','MarkerSize',1)
plot(BiTeSbTeSLOxxx, y,'.r','MarkerSize',1)
hold off
xlabel('Delay in ps','Fontsize',13)
ylabel('\DeltaR / R_0','Fontsize',13)
hleg1=legend('Gemessen','Fit-Funktion, s.Text');
hLine = findobj(hleg1,'type','line','Marker','.');
set(hLine,'MarkerSize',5);
grid

s4=subplot(2,2,4);
cc = flipud(autumn(wvlsize));
h = zeros(3,wvlsize);
hold on
for i = 1:wvlsize;
h(:,i)=plot(fftxxmatrix(:,1),fftyymatrix(:,i),'color',cc(i,:),'markersize',fftmarkersize,'Marker','.');
set(h(:,i),'Color',cc(i,:));
end
hold off
xlabel('\nu in THz','fontsize',13)
ylabel('Amplitudendichte |Y|','fontsize',13)
axis xy;
axis(s4,[0 5 0 1.1e-5]);
grid
print(fig, '26072012BiTeSbTe(10_50)10000all', '-depsc');
