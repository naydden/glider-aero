% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 6;
Uinf = [1,0,0];

A=linspace(0.75*A_ratio,1.25*A_ratio,21);

MGC=0.5*(cr_W+ct_W);
b_W=A*MGC;
sweep_W=0; dihedral_W=0; twist_W=0; 
x_offset_W=0; z_offset_W=1*MGC;

L=zeros(1,length(A));
Dind=zeros(1,length(A));
M=zeros(1,length(A));
CL=zeros(1,length(A));
CD=zeros(1,length(A));
Cm=zeros(1,length(A));

for i=1:length(A)
    % Wing
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W(i),Nx,...
        Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);

    %Plane incidence
    incidence=alpha;
    [Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);

    %Symetric plane
    [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);

    %Ground-effect assembly
    [Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

    % COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    deltaY = b_W(i)/(2*Ny);

    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');

    L(i) = lift(dLw,dLh,dLv);
    M(i) = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind(i) = drag(dDw,dDh,dDv); 
    CDparw = cdragpar(dLw,deltaY(1),Ny,cr_W,ct_W,b_W(i),rho,Uinf,'ala');
    CDpar = [CDparw 0 0];

    [CL(i), CD(i), Cm(i)] = Coeff(cr_W,ct_W,b_W(i),Uinf,rho,L(i),Dind(i),CDpar,M(i));

end

fprintf('Wing case + ground: L= %f D= %f M= %f\n CL= %f CD=%f Cm=%f \n',L(1),Dind(1),M(1),CL(1),CD(1),Cm(1))

% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
surf(Coord(:,1:2*(Ny+1),1),Coord(:,1:2*(Ny+1),2),Coord(:,1:2*(Ny+1),3));
hold on;
surf(Coord(:,2*(Ny+1)+1:4*(Ny+1),1),Coord(:,2*(Ny+1)+1:4*(Ny+1),2),Coord(:,2*(Ny+1)+1:4*(Ny+1),3));
axis equal;

figure(2);
plot(A,CL);
xlabel('Allargament');
ylabel('CL');
grid on;
% createfigure(A, CL, 'A', 'Allargament', 'CL', 'CL en funció del allargament','CL_A');

figure(3);
plot(A,CD);
xlabel('Allargament');
ylabel('CD');
grid on;
% createfigure(A, CD, 'A', 'Allargament', 'CD', 'CD en funció del allargament','CD_A');

figure(4);
plot(A,Cm);
xlabel('Allargament');
ylabel('Cm');
grid on;
% createfigure(A, Cm, 'A', 'Allargament', 'Cm', 'Cm en funció del allargament','Cm_A');