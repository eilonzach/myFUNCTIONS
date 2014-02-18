function [ out_vect ] = midpts( in_vect )
% [ out_vect ] = midpts( in_vect )
%   This function takes a vector of length N and returns a vector of length
%   N-1 consisting of the linear-interpolated midpoints of the N values in
%   the input vector
%
% at some point expand this function to do 2,3 DIM matrices


    out_vect = 0.5*(in_vect(1:end-1)+in_vect(2:end));


end

