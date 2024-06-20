function a_normed = normalise(a)

a_normed = (a - min(a))./(max(a) - min(a));

end

