% Function that returns the momement about the origin of coordinates
% of the wing given the lift of each element and its position.
function M = moment(dLw,dLh,dLv,Nx,Ny,Xp,cas)
    Mw = 0;
    Mh = 0;
    Mv = 0;
    switch cas
        case 'ala'
            Xp_w = rearrange_wing(Nx,Ny,Xp,'wing');
            for j=1:2*Ny
                for i=1:Nx
                  Mw = Mw + Xp_w(i,j)*dLw(i,j);  
                end
            end
        case 'ala+htp'
            Xp_w = rearrange_wing(Nx,Ny,Xp,'wing');
            Xp_h = rearrange_wing(Nx,Ny,Xp,'htp');
            for j=1:2*Ny
                for i=1:Nx             
                  Mw = Mw + Xp_w(i,j)*dLw(i,j);
                  Mh = Mh + Xp_h(i,j)*dLh(i,j);                    
                end
            end            
        case 'ala+htp+vtp'
            Xp_w = rearrange_wing(Nx,Ny,Xp,'wing');
            Xp_h = rearrange_wing(Nx,Ny,Xp,'htp');
            Xp_v = rearrange_wing(Nx,Ny,Xp,'vtp');
            for j=1:Ny
                for i=1:Nx
                      Mw = Mw + Xp_w(i,j)*dLw(i,j);
                      Mh = Mh + Xp_h(i,j)*dLh(i,j);
                      if(j>Ny)
                          continue                         
                      else
                        Mv = Mv + Xp_v(i,j)*dLv(i,j); 
                      end
                end
            end            
    end
    M = Mw + Mh + Mv;
    M = -M;
end
