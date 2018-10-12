function Par = Hyper(XY)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit with "hyperaccuracy" (with zero essential bias)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for speed, not for stability
%
%--------------------------------------------------------------------------

X = XY(:,1);
Y = XY(:,2);
Z = X.*X + Y.*Y;
ZXY1 = [Z X Y ones(length(Z),1)];
M = ZXY1'*ZXY1;
S = mean(ZXY1);
N = [8*S(1) 4*S(2) 4*S(3) 2; 4*S(2) 1 0 0; 4*S(3) 0 1 0; 2 0 0 0];
NM = inv(N)*M;
[E,D] = eig(NM);
[Dsort,ID] = sort(diag(D));
if (Dsort(1)>0) 
    disp('Error in Hyper: the smallest e-value is positive...')
end
if (Dsort(2)<0) 
    disp('Error in Hyper: the second smallest e-value is negative...')
end
A = E(:,ID(2));

Par = zeros(1,3);
Par(1:2) = -(A(2:3))'/A(1)/2;
Par(3) = sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2;

end   %  Hyper