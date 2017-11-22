
% Function that returns the drag of the wing given the drag of each
% element

function Dind = drag(dDind,Nx,Ny)
    Dind = 0;
    for j=1:Ny
        for i=2:Nx
            Dind = Dind + dDind(i,j);
        end
    end
end