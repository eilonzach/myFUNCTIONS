function [ cc,CC6 ] = Etensor_hex( A,C,F,L,N,ax )
%[ cc,CC6 ] = elastic_constants( A,C,F,L,N )
%   function to create the elastic tensor from the 5 Love parameters
% ax denotes which is the symmetry axis - by default is 3, but seismos
% like 1

if ax == 3
CC6 = [   A   A-2*N    F    0  0  0
         A-2*N   A     F    0  0  0
           F     F     C    0  0  0
           0     0     0    L  0  0
           0     0     0    0  L  0
           0     0     0    0  0  N ];
elseif ax == 1
 CC6 = [   C     F     F    0  0  0
           F     A   A-2*N  0  0  0
           F   A-2*N   A    0  0  0
           0     0     0    N  0  0
           0     0     0    0  L  0
           0     0     0    0  0  L ];   
elseif ax == 11
 CC6 = [   C     F   A-2*N  0  0  0
           F     C   A-2*N  0  0  0
         A-2*N A-2*N   A    0  0  0
           0     0     0    N  0  0
           0     0     0    0  N  0
           0     0     0    0  0  L ];   
end
           
           
cc=zeros(3,3,3,3);  % full elasticity
for i=1:3
    for j=1:3
        cc(i,i,j,j)=CC6(i,j);
    end
end
for i=1:2
    for j=(i+1):3
        k=9-j-i;
        c=CC6(k,k);
        cc(i,j,i,j)=c;
        cc(i,j,j,i)=c;
        cc(j,i,j,i)=c;
        cc(j,i,i,j)=c;     
    end
end

end

