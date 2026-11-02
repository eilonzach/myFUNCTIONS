function ax = subplotij(figh,Ny,Nx,i,j)
% ax = subplotij([figh],Ny,Nx,i,j)
% OR
% ax = subplotij([figh],Ny,Nx,i,AX) 
% 
% function to place a subplot in a specific (i,j) location in a grid (saves
% having to work out how to go from row/col to subplot position. 
% If the first argument is a figure handle, will insert into that figure;
% If the last arg [AX or j] is an axes array with Ny rows and Nx columns,
% for instance made with subplot_custom_make, then the output of this 
% function is the axes handle for the i-th subplot, going across rows
% first.
% 
% see also: subplot_custom_make
if nargin < 5
    j = 1;
end

if ishandle(figh) % prexisting figure handle
else % shuffle all the inputs back
    j=i;
    i=Nx;
    Nx=Ny;
    Ny=figh;
    figh = gcf;
end


if ishandle(j)
    AX = j;
    if i > numel(AX), error('Not this many elements in AX'); end
    ii = mod(i-1,Nx)+1;
    jj = ceil(i/Nx);
    ax = AX(jj,ii);
else
    ax = subplot(figh,Ny,Nx,j + (i-1)*Nx);
end

    


end