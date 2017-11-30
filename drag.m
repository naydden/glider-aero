
%% Function that returns the drag of the wing given the drag of each
% element
function Dind = drag(dDw,dDh,dDv)
    Dind= sum(sum(dDw))+sum(sum(dDh))+sum(sum(dDv));
end