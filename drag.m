
% Function that returns the drag of the wing given the drag of each
% element

function D = drag(dD,Nx,Ny)
    D = 0;
    for j=1:Ny
        for i=2:Nx
            D = D + dD(i,j);
        end
    end
end