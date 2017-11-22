% Function that returns the momement about the origin of coordinates
% of the wing given the lift of each element and its position.
function M = moment(dL,Nx,Ny,Xp)
    M = 0;
    for j=1:Ny
        for i=1:Nx
                 
          Xpm = (Xp(i,j)+Xp(i,j+1))/2;
          M = M + Xpm*dL(i,j);
          
        end
    end
    M = -M;
end