clear all;
printf( '=======\n' )

R_obj = expSo3( [0.7;0.3;0.5] );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ]

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 8; 5; 1], ... % necessary to avoid falling in, what seems to be, a local minimum
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;
ys = Ys( 1:2,: ) ./ Ys( 3,: );

n_dims = size( ys, 1 );
n_pts  = size( ys, 2 );


[ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();
J = zeros( n_dims*n_pts, 6 );

T = eye(4);

n_iter = 15
step   = 0.7;
for iter=1:n_iter

  Txs = T*Xs;
  txs = Txs(1:2,:) ./ Txs(3,:);
  e   = ( txs - ys )(:);
  norms(iter) = norm(e);

  for ii=1:n_pts

    d_proj = 1/Txs(3,ii) * [ [ 1, 0, -txs(1,ii), 0 ]; ...
                             [ 0, 1, -txs(2,ii), 0 ] ];

    i_one = ii*2 - 1;
    i_snd = ii*2;

    J( i_one:i_snd, 1 ) = ( d_proj * G_yz * Txs(:,ii) )(:); % R_yz
    J( i_one:i_snd, 2 ) = ( d_proj * G_zx * Txs(:,ii) )(:); % R_zx
    J( i_one:i_snd, 3 ) = ( d_proj * G_xy * Txs(:,ii) )(:); % R_xy
    J( i_one:i_snd, 4 ) = ( d_proj * G_x  * Txs(:,ii) )(:); % t_x
    J( i_one:i_snd, 5 ) = ( d_proj * G_y  * Txs(:,ii) )(:); % t_y
    J( i_one:i_snd, 6 ) = ( d_proj * G_z  * Txs(:,ii) )(:); % t_z

  end

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ ); % or Levenberg-Marquardt
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step * delta;

  T = expSe3(delta) * T;

end

printf( '----------\n' )
Txs = T*Xs;
txs = Txs( 1:2,: ) ./ Txs( 3,: );
e   = ( txs - ys )( : );
norms(iter+1) = norm(e);
T*inv( T_obj )
plot( norms )
