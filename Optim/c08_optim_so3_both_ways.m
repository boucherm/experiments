clear all;
printf('=======\n')

R_obj = [ 0 -1 0;
          1  0 0;
          0  0 1 ];

Xs    = [ [ 1, 0 ];
          [ 0, 0 ];
          [ 1, 1 ] ];
Ys    = R_obj * Xs;
n_pts = size( Xs, 2 );

[ G_yz, G_zx, G_xy ] = so3Generators();

R = eye(3);

n_iter = 5;
step   = 1;
for iter=1:n_iter

  e = [];
  J = [];
  for index=1:n_pts
    X = Xs( :, index );
    Y = Ys( :, index );
    e = [ e; R*X - Y ];

    j( 1:3, 1 ) = G_yz*R*X;
    j( 1:3, 2 ) = G_zx*R*X;
    j( 1:3, 3 ) = G_xy*R*X;
    J = [ J; j ];
  end
  for index=1:n_pts
    X = Xs( :, index );
    Y = Ys( :, index );
    e = [ e; X - R'*Y ];

    j( 1:3, 1 ) = R'*G_yz*Y;
    j( 1:3, 2 ) = R'*G_zx*Y;
    j( 1:3, 3 ) = R'*G_xy*Y;
    J = [ J; j ];
  end
  norms(iter) = norm(e);

  JtJ  = J'*J;
  iJtJ = pinv(JtJ);
  Jte  = J'*e;

  delta = iJtJ * -Jte;
  delta = step*delta;
  R     = expSo3( delta ) * R;

end

printf( '--- \n' )
e = [];
for index=1:size( Xs, 2 )
  X = Xs( :, index );
  Y = Ys( :, index );
  e = [ e; R*X - Y ];
end
norms(iter+1) = norm(e);
R'*R_obj

plot( norms )
