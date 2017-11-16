function cd = cdrag(naca,cl)
if naca == '2412'
    cd = 0.0063 - 0.0033*cl + 0.0067*cl^2;
elseif naca == '0009'
    cd = 0.0055 + 0.0045*cl^2;
else
    error('Input must be 2412 or 0009')
end