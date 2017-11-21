% Function that returns the lift of the wing given the lift of each
% element.
function L = lift(dL,Nx,Ny)
    L = 0;
    for j=1:Ny
        for i=2:Nx
            L = L + dL(i,j);
        end
    end
end