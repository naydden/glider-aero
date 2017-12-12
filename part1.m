twist_angle = 0:-1:-8;
miau = size(twist_angle,2);
alpha_angle = zeros(1,miau);
CD0 = zeros(1,miau);

for i = 1:miau
    [alpha_angle(i), CD0(i)] = ZLangle(cr_W,ct_W,b_W,Nx,Ny,m_W,p_W,sweep_W,dihedral_W,twist_angle(i),x_offset_W,z_offset_W,rho);
end

figure(1);
plot(twist_angle, alpha_angle);
xlabel('Twist angle (�)')
ylabel('\alpha_{ZL} (�)')
grid on;

figure(2);
plot(twist_angle, CD0);
xlabel('Twist angle (�)')
ylabel('C_{D}')
axis([-8 0 0 0.010])
grid on;