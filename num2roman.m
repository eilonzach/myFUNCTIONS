function R = num2roman(x)
    xin = x; % save
    d = [ 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
    c =  {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
    R = [];
    for i = 1:numel(d)
        if x >= d(i)
            z = floor(x/d(i));
            R = [R, repmat(c{i}, 1,z)];
            x = x - z*d(i); % keep only remainder
        end
    end
end