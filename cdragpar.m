function CDpar = cdragpar(dL,deltaY,Ny,cr,ct,b,rho,Uinf,cas)
    S = b*(cr+ct)/2;
    y=linspace(0,b/2,Ny+1); %Y Coordinate of the points of the elemnts
    c=cr-(cr-ct)*y/(b/2); %Computation of the chord of the wing segment
    q = 0.5*rho*norm(Uinf)^2;    
    switch cas
        case 'ala'
            CDpar = 0;
            c_new = zeros(1,Ny);
            for j=1:(Ny)
                c_new(j) = mean([c(j) c(j+1)]);
            end
            % chord of the whole wing 
            c = [fliplr(c_new) c_new];
            for j=1:2*(Ny)
                l = sum(dL(:,j));
                Cl = l/(q*c(j));
                Cd = cdrag('2412',Cl);
                CDpar = CDpar + Cd*c(j)*deltaY/S;
            end
        case 'htp'
            CDpar = 0;
            c_new = zeros(1,Ny);
            for j=1:(Ny)
                c_new(j) = mean([c(j) c(j+1)]);
            end
            % chord of the whole wing 
            c = [fliplr(c_new) c_new];        
            for j=1:2*(Ny)
                l = sum(dL(:,j));
                Cl = l/(q*c(j));
                Cd = cdrag('0009',Cl);
                CDpar = CDpar + Cd*c(j)*deltaY/S;
            end
        case 'vtp'
            CDpar = 0;
            c_new = zeros(1,Ny);
            for j=1:(Ny)
                c_new(j) = mean([c(j) c(j+1)]);
            end            
            for j=1:(Ny)
                l = sum(dL(:,j));
                Cl = l/(q*c_new(j));
                Cd = cdrag('0009',Cl);
                CDpar = CDpar + Cd*c_new(j)*deltaY/S;
            end            
    end   
end