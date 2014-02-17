function [ out_vect ] = midpts( in_vect )
% [ out_vect ] = midpts( in_vect )
%   This function takes a vector of length N and returns a vector of length
%   N-1 consisting of the linear-interpolated midpoints of the N values in
%   the input vector
%
% at some point expand this function to do 2,3 DIM matrices

N = length(in_vect);

if size(in_vect,1) > size(in_vect,2)
out_vect = zeros(N-1,1);
else out_vect = zeros(1,N-1);
end

for ii = 1:N-1
    out_vect(ii) = mean([in_vect(ii),in_vect(ii+1)]);
end


end

