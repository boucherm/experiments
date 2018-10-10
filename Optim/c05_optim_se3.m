clear all;
printf('=======\n')

R_obj = expSo3( [0.7;0.3;0.5] );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ]

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;

n_pts = size( Xs, 2 );

[ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();

T = eye( 4 );

n_iter = 5
step   = 1;
for iter=1:n_iter

  e = [];
  J = [];
  for index=1:n_pts
    x = Xs( :, index );
    y = Ys( :, index );

    e = [ e; T*x - y ];

    j( 1:4, 1 ) = G_yz * T * x; % R_yz
    j( 1:4, 2 ) = G_zx * T * x; % R_zx
    j( 1:4, 3 ) = G_xy * T * x; % R_xy
    j( 1:4, 4 ) = G_x  * T * x; % t_x
    j( 1:4, 5 ) = G_y  * T * x; % t_y
    j( 1:4, 6 ) = G_z  * T * x; % t_z

    J = [ J; j ];
  end
  norms(iter) = norm(e);

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ ); % or Levenberg-Marquardt
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step * delta;

  T = expSe3(delta)*T;

end

printf( '--- \n' )
e = [];
for index=1:n_pts
  x = Xs( :, index );
  y = Ys( :, index );

  e = [ e; T*x - y ];
end
norms(iter+1) = norm(e);

T
inv(T) * T_obj
plot( norms )
