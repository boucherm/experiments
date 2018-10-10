% Perform matrix exponential of an so(3) algebra vector
% IN
%  w: vector in so(3) lie algebra
% OUT
%  R: the SO(3) group resulting matrix

function R = expSo3( w )

  wx     = hat(w);
  theta  = norm(w);
  theta2 = theta * theta;
  if ( abs( theta ) > 0.0000000001 )
    a = sin(theta) / theta;
    b = ( 1 - cos(theta) ) / theta2;
  else
    a = 1;
    b = 0.5;
  end

  R = eye(3) + a*wx + b*wx*wx;

end
