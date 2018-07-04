function [ y ] = window_1(width, EPoW, t0, t)
% EPoW = Edge: Percent of Window
T=EPoW*width.*2;
phi=2*pi./T.*(t0-width./2);
phi2=2*pi./T.*(t0+width./2);
y=zeros(size(t));

a1=t0-(width./2-T./2);
a2=t0+(width./2-T./2);
b1=t0-width./2;
b2=t0+width./2;
for i=1:length(t)
    if (t(i)>=a1) && (t(i)<=a2)
        y(i)=1;
    elseif (t(i)>b1) && (t(i)<a1)
        y(i)=(1-cos(2.*pi./T.*t(i)-phi))./2;
    elseif (t(i)>a1) && (t(i)<b2)
        y(i)=(1-cos(2.*pi./T.*t(i)-phi2))./2;
    else
        y(i)=0;
    end
end

