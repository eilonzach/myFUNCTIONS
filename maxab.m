function [ B ] = maxab( A )
% [ B ] = maxab( A )
%   This function returns the largest value in data series A (or in each
%   column, if A is a matrix. This will be the largest absolute value,
%   irrespective if that is negative or positive

mabA = max(abs(A));
abA = abs(A);
B = zeros(1,size(A,2));
for ii = 1:size(A,2)
    B(ii) = unique(A(abA(:,ii)==mabA(ii),ii));
end


end

