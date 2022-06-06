function [y,I,J] = uniqcount(x)
% UNIQCOUNT Count the number of occurrences of each unique row in matrix
%
%   y = UNIQCOUNT(x) returns a matrix y with one row for each unique row 
%   in x and an additional column to the right giving the count (i.e., the
%   number of occurrences) for each of the unique rows.
%
%   [y,I,J] = UNIQCOUNT(x) in addition returns index vectors I,J such that 
%   y(I,:) == x(J,1:(end-1)) 


    [dum, I, J] = unique(x,'rows');
    y = accumarray(J,1);
    y = [dum y];



  