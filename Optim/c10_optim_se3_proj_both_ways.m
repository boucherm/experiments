clear all;
printf( '=======\n' )

R_obj = expm( hat( [0.7;0.3;0.5] ) );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ];

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       %[0; 8; 5; 1], ... % not necessary anymore ( see c07_optim_se3_proj.m )
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;
xs = Xs(1:2,:) ./ Xs(3,:);
ys = Ys(1:2,:) ./ Ys(3,:);

n_pts  = size( ys, 2 );
n_dims = size( ys, 1 );

[ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();
J_fwd = zeros( n_dims*n_pts, 6 );
e_fwd = zeros( n_dims*n_pts, 1 );
J_bwd = zeros( n_dims*n_pts, 6 );
e_bwd = zeros( n_dims*n_pts, 1 );

T = eye(4);

n_iter = 15;
step   = 0.7;
for iter=1:n_iter

  % Build "forward" error and jacobian
  Txs   = T*Xs;
  txs   = Txs( 1:2,: ) ./ Txs( 3,: );
  e_fwd = ( txs - ys )( : );

  for ii=1:n_pts

    d_proj = 1/Txs(3,ii) * [ [ 1, 0, -txs(1,ii), 0 ]; ...
                             [ 0, 1, -txs(2,ii), 0 ] ];

    i_one = ii*2 - 1;
    i_snd = ii*2;

    J_fwd( i_one:i_snd, 1 ) = ( d_proj * G_yz * Txs(:,ii) )(:); % R_yz
    J_fwd( i_one:i_snd, 2 ) = ( d_proj * G_zx * Txs(:,ii) )(:); % R_zx
    J_fwd( i_one:i_snd, 3 ) = ( d_proj * G_xy * Txs(:,ii) )(:); % R_xy
    J_fwd( i_one:i_snd, 4 ) = ( d_proj * G_x  * Txs(:,ii) )(:); % t_x
    J_fwd( i_one:i_snd, 5 ) = ( d_proj * G_y  * Txs(:,ii) )(:); % t_y
    J_fwd( i_one:i_snd, 6 ) = ( d_proj * G_z  * Txs(:,ii) )(:); % t_z

  end

  % Build "backward" error and jacobian
  Tinv = inv(T);

  Tiys  = Tinv*Ys;
  tiys  = Tiys(1:2,:) ./ Tiys(3,:);
  e_bwd = ( tiys - xs )(:);

  for ii=1:n_pts

    d_proj = 1/Tiys(3,ii) * [ [ 1, 0, -tiys(1,ii), 0 ]; ...
                              [ 0, 1, -tiys(2,ii), 0 ] ];
    i_one = ii*2 - 1;
    i_snd = ii*2;

    J_bwd( i_one:i_snd, 1 ) = ( d_proj * Tinv * -G_yz * Ys(:,ii) )(:); % R_yz
    J_bwd( i_one:i_snd, 2 ) = ( d_proj * Tinv * -G_zx * Ys(:,ii) )(:); % R_zx
    J_bwd( i_one:i_snd, 3 ) = ( d_proj * Tinv * -G_xy * Ys(:,ii) )(:); % R_xy
    J_bwd( i_one:i_snd, 4 ) = ( d_proj * Tinv * -G_x  * Ys(:,ii) )(:); % t_x
    J_bwd( i_one:i_snd, 5 ) = ( d_proj * Tinv * -G_y  * Ys(:,ii) )(:); % t_y
    J_bwd( i_one:i_snd, 6 ) = ( d_proj * Tinv * -G_z  * Ys(:,ii) )(:); % t_z

  end

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
Txs = T*Xs;
txs = Txs(1:2,:) ./ Txs(3,:);
e   = ( txs - ys )(:);
norms(iter+1) = norm(e);
T*inv( T_obj )
plot( norms )
