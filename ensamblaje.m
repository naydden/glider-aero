function A = ensamblaje(B,C)
if ndims(B)==3
    A = [ B , C ];
else
    A = [ B ; C ];
end