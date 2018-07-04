function [autocorrtime, autocor] = prep(datatemp, zeroref_start, zeroref_end)
time = datatemp(:,1);
data = datatemp(:,2);

% Daten nur von timemin bis timemax

[~,I1]=min(time);
[~,I2]=max(time);
if I2>I1
    a = 1;
else
    a=-1;
end

autocorrtime = time(I1:a:I2);
autocorr = data(I1:a:I2);

% Zeitachse �quidistant
% doppelte zeitwerte rauskicken
[xxx,yyy] = consolidator(autocorrtime, autocorr, 'mean');
delta = abs(max(xxx)-min(xxx))./length(xxx);
measuretime = (min(xxx):delta:max(xxx))';
%reflektivit�tswerte interpolieren
autocorr1 = interp1(xxx,yyy,measuretime);

% Ref-Achse eichen (plot(measureref1))
meanref = mean(measureref1(zeroref_start:zeroref_end));
measureref = measureref1 - meanref;