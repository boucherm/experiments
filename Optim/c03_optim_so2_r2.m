clear all;
printf('=======\n')

R_obj = [ 0 -1;
          1  0 ]
t_obj = [ 1;
          1 ]

x1 = [ 1;
       0 ]
y1 = R_obj*x1 + t_obj

x2 = [ 0;
       1 ]
y2 = R_obj*x2 + t_obj

G = [ 0 -1;
      1  0 ];

R = eye(2);
t = [ 0; 0 ];

n_iter = 100
step   = 0.1;

for iter=1:n_iter

  e( 1:2, 1 ) = R*x1 + t - y1;
  e( 3:4, 1 ) = R*x2 + t - y2;
  norms(iter) = norm(e);

  J( 1:2, 1 ) = G*R*x1;
  J( 1:2, 2 ) = [ 1; 0 ];
  J( 1:2, 3 ) = [ 0; 1 ];

  J( 3:4, 1 ) = G*R*x2;
  J( 3:4, 2 ) = [ 1; 0 ];
  J( 3:4, 3 ) = [ 0; 1 ];

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ );
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step*delta;

  R    = expm( delta(1) * G ) * R;
  t(1) = t(1) + delta(2);
  t(2) = t(2) + delta(3);

end

printf( '--- \n' )
e( 1:2, 1 ) = R*x1 + t - y1;
e( 3:4, 1 ) = R*x2 + t - y2;
norms(iter+1) = norm(e);

R
t
plot( norms )
