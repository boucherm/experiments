% Perform matrix exponential of an se(3) algebra vector
% IN
%  v: vector in se(3) lie algebra, as [ w; u ] for rotational and translational parts
% OUT
%  T: the SE(3) group resulting matrix

function T = expSe3( v )

  w = v( 1:3 );
  u = v( 4:6 );

  wx     = hat(w);
  theta  = norm(w);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  if ( abs(theta) > 0.0000000001 )
    a = (   1   - cos( theta ) ) / ( theta2 );
    b = ( theta - sin( theta ) ) / ( theta3 );
  else
    a = (1/2) * ( 0.5 - (theta2/12)*( 1 - (theta2/30)*(1-theta2/56) ) );
    b = (1/6) * (  1  - (theta2/20)*( 1 - (theta2/42)*(1-theta2/72) ) );
  end

  V = eye(3) + a*wx + b*wx*wx;
  R = expSo3(w);
  T = [ [ R,  V*u ]; [ 0 0 0 1 ] ]; % TODO define and use expSo3

end
