% Function that returns the lift of the wing given the lift of each
% element.
function L = lift(Nx,Ny,dLw,dLh,dLv)
    L = 0;
    for j=1:Ny
        for i=2:Nx
             L = L + dLw(i,j)+dLh(i,j)+dLv(i,j);
        end
    end
end    
    