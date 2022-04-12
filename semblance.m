function s = semblance(mod1,mod2, ifnorm, thresh)
%s = semblance(mod1,mod2, ifnorm, thresh)
%
%  Calcualte semblance - agreement between two models after Zelt 1998
% 
%  We calculate the mean semblance at each node by summing over it and its
%  neighbours' coherence between input and output, over some radius - make
%  sure the mod1 and mod2 are vectors of points only within some radius of
%  each other
% 
% Z. Eilon 03/2016

if nargin <3
    ifnorm = 0;
end


if ifnorm
   mod1(mod1>thresh) = thresh;
   mod1(mod1<-thresh) = -thresh;
   mod1(abs(mod1)<thresh) = 0;

   mod2(mod2>thresh) = thresh;
   mod2(mod2<-thresh) = -thresh;
   mod2(abs(mod2)<thresh) = 0;
end 


top = nansum( (mod1 + mod2).^2 );
bot = 2*nansum( mod1.^2 + mod2.^2 );
s = top./bot;


end