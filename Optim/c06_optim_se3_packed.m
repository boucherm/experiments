clear all;
printf('=======\n')

R_obj = expSo3( [0.7;0.3;0.5] );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ]

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;

n_dims = size( Xs, 1 );
n_pts  = size( Xs, 2 );

[ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();
J = zeros( n_dims*n_pts, 6 );

T = eye(4)

n_iter = 5;
step   = 1;
for iter=1:n_iter

  e = ( T*Xs - Ys )(:);
  norms(iter) = norm(e);

  Txs = T*Xs;

  J(:,1) = ( G_yz * Txs )(:); % R_yz
  J(:,2) = ( G_zx * Txs )(:); % R_zx
  J(:,3) = ( G_xy * Txs )(:); % R_xy
  J(:,4) = ( G_x  * Txs )(:); % t_x
  J(:,5) = ( G_y  * Txs )(:); % t_y
  J(:,6) = ( G_z  * Txs )(:); % t_z

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ ); % or Levenberg-Marquardt
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step * delta;

  T = expSe3(delta) * T;

end

printf( '----------\n' )
e = ( T*Xs - Ys )(:);
norms(iter+1) = norm(e);
T*inv( T_obj )
plot( norms )
