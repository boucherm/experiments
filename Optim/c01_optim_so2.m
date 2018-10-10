clear all;
printf('=======\n')

R_obj = [ 0 -1; 1 0 ];

x = [1; 0];
y = R_obj * x;

G = [0 -1; 1 0];

R = eye(2);

n_iter = 10
step   = 1;
for iter=1:n_iter

  e = R*x - y;
  norms(iter) = norm(e);

  J = G*R*x; % Update on the left
  %J = R*G*x; % Update on the right

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ );
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step*delta;

  R = expm(delta*G) * R; % Update on the left
  %R = R * expm(delta*G); % Update on the right

end

printf('---\n')
e = R*x - y;
norms(iter+1) = norm(e);

R'*R_obj
plot( norms )
