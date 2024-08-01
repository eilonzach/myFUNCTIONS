function [ height_to_width ] = plot_size_ratio( lonlims,latlims )
% [ height_to_width ] = plot_size_ratio( lonlims,latlims )
%   Function to calculate lat distance vs lon distance to have right
%   dimensions on plots To have correct dimensions in map, the y dimension
%   [ax.Position(4)] of the plot (in pixels) should be equal to the x
%   dimension [ax.Position(3)] MULTIPLIED by the height_to_width value.

height = distance_km(latlims(1),lonlims(1),latlims(2),lonlims(1));
width  = distance_km(latlims(1),lonlims(1),latlims(1),lonlims(2));

height_to_width = height/width;


end