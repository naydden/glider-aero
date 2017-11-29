% Function that returns the lift of the wing given the lift of each
% element.
function L = lift(dLw,dLh,dLv)
    L = sum(sum(dLw))+sum(sum(dLh))+sum(sum(dLv));
end    
    