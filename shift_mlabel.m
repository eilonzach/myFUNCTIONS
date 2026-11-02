function ax =shift_mlabel(ax,shiftfactor,Nstr,Estr,Sstr,Wstr)
% ax = shift_mlabel(ax,shiftfactor,Nstr,Estr,Sstr,Wstr)
% 
% function to shift axes labels off m_axis objects. Will always shift away
% from the center inputs are the axes handle contianing the m_grid, the
% factor by which you want to shift the labels,and then vectors (doubles)
% with the values of any ticks you want to shift from the N, E, S, and W
% sides respectively.

    hhh = get(ax,'Children');
    for ih = 1:length(hhh)
        if ~contains(hhh(ih).Tag,'m_grid'), continue; end
        if ~contains(hhh(ih).Tag,'ticklabel'), continue; end
        for ihn = 1:length(Nstr)
            if contains(hhh(ih).String,'N')
                if contains(hhh(ih).String,num2str(Nstr(ihn)))
                    hhh(ih).Position(1) = hhh(ih).Position(1)*shiftfactor;
                end
            end
        end
        for ihn = 1:length(Sstr)
            if contains(hhh(ih).String,'S')
                if contains(hhh(ih).String,num2str(Sstr(ihn)))
                    hhh(ih).Position(1) = hhh(ih).Position(1)*shiftfactor;
                end
            end
        end
        for ihn = 1:length(Estr)
            if contains(hhh(ih).String,'E')
                if contains(hhh(ih).String,num2str(Estr(ihn)))
                    hhh(ih).Position(2) = hhh(ih).Position(2)*shiftfactor;
                end
            end
        end
        for ihn = 1:length(Wstr)
            if contains(hhh(ih).String,'W')
                if contains(hhh(ih).String,num2str(Wstr(ihn)))
                    hhh(ih).Position(2) = hhh(ih).Position(2)*shiftfactor;
                end
            end
        end
    end


end