function Gamma = rearrange_wing(Nx,Ny,Gamma,cas)
    switch cas
        case 'wing'
            Gamma_e = zeros(Nx,Ny);
            Gamma_d = zeros(Nx,Ny);  
            k = 1;
            % semiala esquerra            
            for i=1:Nx
                for j=1:Ny
                    % change from vector to matrix
                    k = k + 1;
                    Gamma_e(i,j) = Gamma((i-1)*Ny+j);
                end
            end
            % semiala dreta
            Ny_aux = k+Ny-1;
            for i=1:Nx
                l = 0;
                for j=k:Ny_aux
                    l = l + 1;
                    % change from vector to matrix
                    Gamma_d(i,l) = Gamma((i-1)*Ny+j);
                end
            end
            Gamma = [Gamma_e Gamma_d];
        case 'vtp'
            Gamma_new = zeros(Nx,Ny);
            start = Nx*Ny*4;
            for i=1:Nx
                l = 0;
                for j=(start+1):(start+Ny)
                    l = l + 1;
                    % change from vector to matrix
                    Gamma_new(i,l) = Gamma((i-1)*Ny+j);
                end
            end
            Gamma = Gamma_new;
        case 'htp'
            Gamma_e = zeros(Nx,Ny);
            Gamma_d = zeros(Nx,Ny);
            start = Nx*2*Ny;
            k = 1;
            % semiala esquerra            
            for i=1:Nx
                l = 0;
                for j=(start+1):(start+Ny)
                    l = l + 1;
                    % change from vector to matrix
                    k = k + 1;
                    Gamma_e(i,l) = Gamma((i-1)*Ny+j);
                end
            end
            % semiala dreta
            start = start + k -1;
            for i=1:Nx
                l = 0;
                for j=(start+1):(start+Ny)
                    l = l + 1;
                    % change from vector to matrix
                    Gamma_d(i,l) = Gamma((i-1)*Ny+j);
                end
            end
            Gamma = [Gamma_e Gamma_d];         
    end
end            