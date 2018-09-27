clear all;
printf('=======\n')

R_obj = expm( hat( [0.7;0.3;0.5] ) );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ]

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;

T = eye( 4 )

G_yz = [ [ hat([1; 0; 0]),  zeros(3, 1) ]; [ zeros(1, 4)] ];
G_zx = [ [ hat([0; 1; 0]),  zeros(3, 1) ]; [ zeros(1, 4)] ];
G_xy = [ [ hat([0; 0; 1]),  zeros(3, 1) ]; [ zeros(1, 4)] ];

G_x = zeros(4); G_x(1, 4) = 1;
G_y = zeros(4); G_y(2, 4) = 1;
G_z = zeros(4); G_z(3, 4) = 1;

J = zeros( size( Ys, 1 ) * size( Ys, 2 ), 6 );

n_iter = 100
step   = 0.1;

for iter=1:n_iter

  e = ( T*Xs - Ys )(:);
  norms(iter) = norm(e);

  Txs = T*Xs;

  J( :, 1 ) = ( G_yz * Txs )(:); % R_yz
  J( :, 2 ) = ( G_zx * Txs )(:); % R_zx
  J( :, 3 ) = ( G_xy * Txs )(:); % R_xy
  J( :, 4 ) = ( G_x  * Txs )(:); % t_x
  J( :, 5 ) = ( G_y  * Txs )(:); % t_y
  J( :, 6 ) = ( G_z  * Txs )(:); % t_z

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ ); % Levenberg-Marquardt could be used instead
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step * delta;

  w = delta( 1:3 );
  u = delta( 4:6 );

  wx     = hat( w );
  theta  = sqrt( w'*w );
  theta2 = theta*theta;
  theta3 = theta2*theta;
  if ( abs( theta ) > 0.000001 )
    a = (   1   - cos( theta ) ) / ( theta2 );
    b = ( theta - sin( theta ) ) / ( theta3 );
  else
    a = (1/2) * ( 0.5 - (theta2/12)*( 1 - (theta2/30)*( 1 - theta2/56 ) ) );
    b = (1/6) * (  1  - (theta2/20)*( 1 - (theta2/42)*( 1 - theta2/72 ) ) );
  end
  V = eye(3) + a*wx + b*wx*wx;
  incr = [ [ expm( wx ),  V*u ]; [ 0 0 0 1 ] ];

  T = incr*T;

end

printf( '----------\n' )
e = ( T*Xs - Ys )(:);
norms(iter+1) = norm(e);
T
plot( norms )