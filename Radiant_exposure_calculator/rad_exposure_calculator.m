%Power of beams measured at sample in mW
Pspump = 800.0;
Psprobe = 1.874;

%%with laser power in mW
Pout = 1780;
 
%Angle of incidence (angle to perpendicular relative to surface)
pumpalpha= 30;
probealpha= 11;

%FWHM beam waists in micrometer
% (put left and right boundary in the brackets)

pumpsigmax=(111.5+113.5)./2;
pumpsigmay=(78.0+79.0)./2;
probesigmax= (28.3+33.1)./2;
probesigmay=(17.3+20.3)./2;

%Ratio Power at Sample/Pout
pumpratio = Pspump./Pout;
proberatio =Psprobe./Pout;

%repetition Rate of your laser in Hz
omega = 80000000;

%Pulse enegy
Epump = Pspump./omega;
Eprobe = Psprobe./omega;

%Resulting fluence on sample in muJ/cm²
pumpfluence = 1000*(Epump*cos(pumpalpha*pi./180))./(0.25*pi*pumpsigmax*pumpsigmay*1e-8)
probefluence = 1000*(Eprobe*cos(probealpha*pi./180))./(0.25*pi*probesigmax*probesigmay*1e-8)
