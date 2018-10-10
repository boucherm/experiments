clear all;
printf('=======\n')

R_obj = [ 0 -1;
          1  0 ];
t_obj = [ 1;
          0 ];

x = [ 1;
      0 ];
y = R_obj*x + t_obj;

G = [ 0 -1;
      1  0 ];

R = eye(2);
t = [ 0; 0 ];

n_iter = 15
step   = 1;
for iter=1:n_iter

  e = R*x + t - y;
  norms(iter) = norm(e);

  J( :, 1 ) = G*R*x;
  J( :, 2 ) = [ 1; 0 ];

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ );
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step*delta;

  R    = expm( delta(1) * G ) * R;
  t(1) = t(1) + delta(2);

end

printf( '--- \n' )
e = R*x + t - y;
norms(n_iter+1) = norm(e);

R'*R_obj
t-t_obj
plot( norms )
