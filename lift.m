% Function that returns the lift of the wing given the lift of each
% element.
function L = lift(dL,ground)
if ground==0
    L=sum(dL);
else
    L=sum(dL(1:0.5*size(dL)));
end
end