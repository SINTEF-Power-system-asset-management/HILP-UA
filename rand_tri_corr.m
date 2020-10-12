function [y1, y2] = rand_tri_corr(a1,b1,c1,a2,b2,c2,rho)
% RAND_TRI_CORR Sample from correlated bivariate triangular distribution
%
%   [Y1, Y2] = RAND_TRI_CORR(A1,B1,C1,A2,B2,C2,RHO) returns two correlated
%   variables y1 and y2 drawn from triangular distributions with
%   correlation rho. The lower, mode and upper values of the triangular
%   distributions are a, b and c, respectively.

% Draw from bivariate correlated normal distribution
z = mvnrnd([0 0],[1 rho; rho 1],1);

% Transform both variables to uniform distributions
u1 = normcdf(z(1));
u2 = normcdf(z(2));

% Transform both variables from uniform to triangular distributions
y1 = tri_inv(u1,a1,b1,c1);
y2 = tri_inv(u2,a2,b2,c2);

end

function x = tri_inv(u,a,b,c)
    % Inverse transform of the cumulative distribution function, x = F^{-1}(u),
    % for a  triangular distribution with lower, mode and upper values a,
    % b, c, respectively.
    t = (b-a)/(c-a);
    if u <= t
        x = a + sqrt((b-a)*(c-a)*u);
    else
        x = c - sqrt((c-a)*(c-b)*(1-u));
    end
end