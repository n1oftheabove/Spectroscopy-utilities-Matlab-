%needs a time vector and a ref vector, outputs time and ref only between
%shaker maxima
function [xxx, yyy] = timecutter(xx,yy)
[~,I1]=min(xx);
[~,I2]=max(xx);
if I2>I1
    a = 1;
else
    a=-1;
end

xxx = xx(I1:a:I2);
yyy = yy(I1:a:I2);