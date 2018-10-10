clear all;
printf('=======\n')

R_obj = expSo3( [0.7;0.3;0.5] );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ];

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;
n_dims = size( Xs, 1 );
n_pts  = size( Xs, 2 );
n_ways = 2;

[ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();

J_fwd = zeros( n_dims*n_pts, 6 );
e_fwd = zeros( n_dims*n_pts, 1 );
J_bwd = zeros( n_dims*n_pts, 6 );
e_bwd = zeros( n_dims*n_pts, 1 );

T = eye(4)

n_iter = 5;
step   = 1;
for iter=1:n_iter

  % Build "forward" error and jacobian
  Txs = T*Xs;

  e_fwd( :, 1 ) = ( Txs - Ys )(:);

  J_fwd(:,1) = ( G_yz * Txs )(:); % R_yz
  J_fwd(:,2) = ( G_zx * Txs )(:); % R_zx
  J_fwd(:,3) = ( G_xy * Txs )(:); % R_xy
  J_fwd(:,4) = ( G_x  * Txs )(:); % t_x
  J_fwd(:,5) = ( G_y  * Txs )(:); % t_y
  J_fwd(:,6) = ( G_z  * Txs )(:); % t_z

  % Build "backward" error and jacobian
  Tinv = inv(T);

  e_bwd( :, 1 ) = ( Tinv*Ys - Xs )(:);

  J_bwd(:,1) = ( Tinv * -G_yz * Ys )(:); % R_yz
  J_bwd(:,2) = ( Tinv * -G_zx * Ys )(:); % R_zx
  J_bwd(:,3) = ( Tinv * -G_xy * Ys )(:); % R_xy
  J_bwd(:,4) = ( Tinv * -G_x  * Ys )(:); % t_x
  J_bwd(:,5) = ( Tinv * -G_y  * Ys )(:); % t_y
  J_bwd(:,6) = ( Tinv * -G_z  * Ys )(:); % t_z

  % Solve
  e = [ e_fwd; e_bwd ];
  J = [ J_fwd; J_bwd ];

  norms(iter) = norm(e);

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
