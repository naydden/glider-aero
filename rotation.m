function CoordR = rotation(Coord,alpha,cr,x_offset,z_offset)

Coord(:,:,1)=Coord(:,:,1)-(x_offset+0.25*cr);
Coord(:,:,3)=Coord(:,:,3)-z_offset;

for i=1:size(Coord,1)
    for j=1:size(Coord,2)
        CoordR(i,j,1)=Coord(i,j,1)*cosd(alpha)+Coord(i,j,3)*sind(alpha);
        CoordR(i,j,3)=Coord(i,j,3)*cosd(alpha)-Coord(i,j,1)*sind(alpha);
    end
end
CoordR(:,:,1)=CoordR(:,:,1)+(x_offset+0.25*cr);
CoordR(:,:,3)=CoordR(:,:,3)+z_offset;
CoordR(:,:,2)=Coord(:,:,2);
end