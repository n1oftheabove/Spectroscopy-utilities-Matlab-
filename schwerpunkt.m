% calculates the centroid of two vectors with customizable index range (start: schwstart, end: schwend)

% think of the vector xxx as a vector consisting of 1 dimensional space coordinates, and yyy as a vector containing the masses corresponding to that coordinates. Then the centrois is calculated as  sum(xxx_i * yyy_i) / sum(yyy_i)


function [schw] = schwerpunkt(xxx,yyy,schwstart,schwend)
A = zeros(1,schwend-schwstart);
for i = schwstart:schwend;
    A(1,i+1-schwstart)=xxx(i)*yyy(i);
end
B = sum(yyy(schwstart:schwend));
schw = sum(A)./B;
