function y = rand_tri(a,b,c)
% TRI_CORR Sample from triangular distribution
%
%   Y = RAND_TRI(A,B,C) returns a variable Y drawn 
%   from a triangular distribution with lower, mode and upper values A, B 
%   and C, respectively.

t = (b-a)/(c-a);
u = rand(1);
if u <= t
    y = a + sqrt((b-a)*(c-a)*u);
else
    y = c - sqrt((c-a)*(c-b)*(1-u));
end