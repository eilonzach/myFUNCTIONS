function ax = subplotij(Ny,Nx,i,j)
% ax = subplotij(Ny,Nx,i,j)
% function to place a subplot in a specific (i,j) location in a grid (saves
% having to work out how to go from row/col to subplot position
ax = subplot(Ny,Nx,j + (i-1)*Nx);

end