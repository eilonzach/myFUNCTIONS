function x = axlim(ax,val)
%  x = axlim(ax,val)
%  Function to grab an element of the axis [xmin xmax ymin ymax] vector 
%    for a particular set of axes. 
%  I.e. 
%       lims = axis(ax);
%       x = lims(val)
%  This is to allow one-line operations on the axes dimensions.

    if nargin < 2
        val = [];
    end

    if ~ishandle(ax)
        val = ax;
        ax = gca;
    end

    try
        lims = axis(ax);
    catch % note will fail if different limit types
        lims = {xlim(ax),ylim(ax)};
        warning('NOTE x and y have different data types. Storing in separate elements of cell')
    end

    if isempty(val)
        val = 1:length(lims);
    end

    x = lims(val);

end


