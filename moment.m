% Function that returns the momement about the origin of coordinates
% of the wing given the lift of each element and its position.
function M = moment(dL,Vortex,ground)
x=0.5*(Vortex(1,:,1)+Vortex(2,:,1))';
dM=-dL.*x;
if ground==0
    M=sum(dM);
else
    M=sum(dM(1:0.5*size(dM)));
end
end