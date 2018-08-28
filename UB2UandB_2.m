%computes U and B from the factor UB

function [U, B] = UB2UandB_2(UB)
BB = UB'*UB;
b11 = sqrt(BB(1,1));
b12 = BB(1,2)./b11;
b13 = BB(1,3)./b11;
b22 = sqrt(BB(2,2)-b12^2);
b23 = (BB(2,3)-b12*b13)/b22;
b33 = sqrt(BB(3,3)- b13^2 - b23^2);

B = [b11 b12 b13;
    0 b22 b23;
    0 0 b33];

U = (UB)*inv(B);

end
