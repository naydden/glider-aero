
%% Function that returns the drag of the wing given the drag of each
% element

function Dind = drag(dDind,Nx,Ny)

    Dind = 0;
    N = length(dDind); %Indicates the number of columns of the matrix
    
    for i=1:Nx
        for j=1:N
            Dind = Dind + dDind(i,j);
        end
    end
end