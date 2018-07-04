%needs xx and yy values, creates a vector xxx from xx with aequidistant values,
%same range and length, and uses this new vector to interpolate
%corresponding yyy-values

function [xxx,yyy] = aequidisttime(xx,yy)

[xx1,yy1] = consolidator(xx, yy, 'mean');
delta = abs(max(xx1)-min(xx1))./length(xx1);
xxx = (min(xx1):delta:max(xx1))';
%reflektivit�tswerte interpolieren
yyy = interp1(xx1,yy1,xxx);


% [xxx,yyy] = consolidator(timecut, corcut, 'mean');
% delta = abs(max(xxx)-min(xxx))./length(xxx);
% measuretime = (min(xxx):delta:max(xxx))';
% %reflektivit�tswerte interpolieren
% measureref1 = interp1(xxx,yyy,measuretime);