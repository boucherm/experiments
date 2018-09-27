clear all;
printf('=======\n')

R_obj = [ 0 -1 0;
          1  0 0;
          0  0 1 ]

Xs = [ [ 1; 0; 1 ], [ 0; 0; 1 ] ]
Ys = R_obj * Xs

R = eye(3);

n_iter = 100
step   = 0.1;

for iter=1:n_iter

  e = [];
  J = [];
  for index=1:size( Xs, 2 )
    x = Xs( :, index );
    y = Ys( :, index );
    e = [ e; R*x - y ];

    j( 1:3, 1 ) = hat([1; 0; 0])*R*x;
    j( 1:3, 2 ) = hat([0; 1; 0])*R*x;
    j( 1:3, 3 ) = hat([0; 0; 1])*R*x;
    J = [ J; j ];
  end
  norms(iter) = norm(e);

  JtJ  = J'*J;
  iJtJ = pinv(JtJ);
  Jte  = J'*e;

  delta = iJtJ * -Jte;
  delta = step*delta;
  R     = expm( hat(delta) ) * R;
end

printf( '--- \n' )
e = [];
for index=1:size( Xs, 2 )
  x = Xs( :, index );
  y = Ys( :, index );
  e = [ e; R*x - y ];
end
norms(iter+1) = norm(e);

R
plot( norms )
