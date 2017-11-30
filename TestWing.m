function tests = TestWing
tests = functiontests(localfunctions);
end
%% DATA
function testCL(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 3; x_offset_W = 0; z_offset_W = 0; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
        Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [CL, ~, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = CL;
    % XFLR5 solution for this data
    expSolution = 0.4261;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
function testCD(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 3; x_offset_W = 0; z_offset_W = 0; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
        Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [~, CD, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = CD;
    % XFLR5 solution for this data
    expSolution = 0.0059;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
function testCm(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 3; x_offset_W = 0; z_offset_W = 0; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
        Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [~, ~, Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = Cm;
    % XFLR5 solution for this data
    expSolution = -0.1572;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
%%%%%%%%%%%%%%%%%%%%% ground
function testCL_ground(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 0; x_offset_W = 0; z_offset_W=1; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    %Plane incidence
    incidence=3;
    [Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);
    %Symetric plane
    [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);
    %Ground-effect assembly
    [Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [CL, ~, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = CL;
    % XFLR5 solution for this data
    expSolution = 0.4734;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
function testCD_ground(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 0; x_offset_W = 0; z_offset_W=1; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    %Plane incidence
    incidence=3;
    [Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);
    %Symetric plane
    [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);
    %Ground-effect assembly
    [Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [~, CD, ~] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = CD;
    % XFLR5 solution for this data
    expSolution = 0.0041;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end
function testCm_ground(testCase)
    % Profile geometry
    m_W = 0.02; p_W = 0.4;
    % Wing geometry
    cr_W = 1; ct_W = 1; b_W = 10; sweep_W = 0; dihedral_W = 0; twist_W = 0;
    % Air
    alpha = 0; x_offset_W = 0; z_offset_W=1; rho = 1.225;
    Uinf = [1*cosd(alpha),0,1*sind(alpha)];
    % Numerical
    Nx = 5; Ny = 5;
    deltaY = b_W/(2*Ny);
    [Coord,Vortex,ControlP,DragP,Normal] = wing_assembly (cr_W,ct_W,b_W,...
    Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_W,x_offset_W,z_offset_W);
    %Plane incidence
    incidence=3;
    [Coord,Vortex,ControlP,DragP,Normal] = rotation(Coord,Vortex,ControlP,DragP,Normal,incidence,0,cr_W,x_offset_W,z_offset_W);
    %Symetric plane
    [Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr] = mirror (Coord,Vortex,ControlP,DragP,Normal);
    %Ground-effect assembly
    [Coord,Vortex,ControlP,DragP,Normal] = assembly(Coord,Vortex,ControlP,DragP,Normal,Coord_Mirr,Vortex_Mirr,ControlP_Mirr,DragP_Mirr,Normal_Mirr);

    Gamma = circulation(Uinf,Vortex,ControlP,Normal);
    [dLw,dLh,dLv] = delta_lift(Gamma,deltaY,Nx,Ny,rho,Uinf,'ala');
    [dDw,dDh,dDv] = delta_drag(Gamma,Vortex,DragP,deltaY,Nx,Ny,rho,Uinf,'ala');
    L = lift(dLw,dLh,dLv);
    M = moment(dLw,dLh,dLv,Nx,Ny,DragP(:,:,1),'ala');
    Dind = drag(dDw,dDh,dDv);
    [~, ~,Cm] = Coeff(cr_W,ct_W,b_W,Uinf,rho,L,Dind,M);
    % Solution
    actSolution = Cm;
    % XFLR5 solution for this data
    expSolution = -0.1712;
    verifyEqual(testCase,actSolution,expSolution,'AbsTol',0.0001)
end