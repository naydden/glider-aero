% Function that returns the momement about the origin of coordinates
% of the wing given the lift of each element and its position.
function M = moment(dL,Nx,Ny,Xp)
    M = 0;
    for j=1:Ny
        for i=2:Nx
        % Not sure how to treat last element             
            if j==Ny
                Xpm = Xp(i,j);
            else
                Xpm = (Xp(i,j)+Xp(i,j+1))/2;
            end
            M = M + Xpm*dL(i,j);
        end
    end
    M = -M;
end