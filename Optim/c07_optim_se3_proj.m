clear all;
printf( '=======\n' )

R_obj = expm( hat( [0.7;0.3;0.5] ) );
t_obj = [ 0.5; -0.9; 0.1 ];
T_obj = [ [ R_obj,  t_obj]; [ 0 0 0 1 ] ]

Xs = [ [1; 0; 1; 1], ...
       [0; 1; 1; 1], ...
       [0; 8; 5; 1], ... % necessary to avoid falling in, what seems to be, a local minimum
       [0; 0; 1; 1] ];
Ys = T_obj * Xs;
ys = Ys( 1:2,: ) ./ Ys( 3,: );

T = eye( 4 );

G_yz = [ [ hat( [1; 0; 0] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];
G_zx = [ [ hat( [0; 1; 0] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];
G_xy = [ [ hat( [0; 0; 1] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];

G_x = zeros( 4 ); G_x( 1, 4 ) = 1;
G_y = zeros( 4 ); G_y( 2, 4 ) = 1;
G_z = zeros( 4 ); G_z( 3, 4 ) = 1;

J = zeros( size( ys, 1 ) * size( ys, 2 ), 6 );

n_iter = 100
step   = 0.1;
for iter=1:n_iter

  Txs = T*Xs;
  txs = Txs( 1:2,: ) ./ Txs( 3,: );
  e   = ( txs - ys )( : );
  norms( iter ) = norm( e );

  for ii=1:size( Xs, 2 )

    % Projection function derivative, see Eade's thesis
    d_proj = 1/Txs( 3,ii ) * [ [ 1 0 -Txs( 1,ii )/Txs( 3,ii ) 0 ]; ...
                               [ 0 1 -Txs( 2,ii )/Txs( 3,ii ) 0 ] ];

    i_one = ii*2 - 1;
    i_snd = ii*2;

    J( i_one:i_snd, 1 ) = ( d_proj * G_yz * Txs( :, ii ) )( : ); % R_yz
    J( i_one:i_snd, 2 ) = ( d_proj * G_zx * Txs( :, ii ) )( : ); % R_zx
    J( i_one:i_snd, 3 ) = ( d_proj * G_xy * Txs( :, ii ) )( : ); % R_xy
    J( i_one:i_snd, 4 ) = ( d_proj * G_x  * Txs( :, ii ) )( : ); % t_x
    J( i_one:i_snd, 5 ) = ( d_proj * G_y  * Txs( :, ii ) )( : ); % t_y
    J( i_one:i_snd, 6 ) = ( d_proj * G_z  * Txs( :, ii ) )( : ); % t_z

  end

  JtJ   = J'*J;
  iJtJ  = pinv( JtJ ); % Levenberg-Marquardt could be used instead
  Jte   = J'*e;
  delta = iJtJ * -Jte;
  delta = step * delta;

  w = delta( 1:3 );
  u = delta( 4:6 );

  wx     = hat( w );
  theta  = sqrt( w'*w );
  theta2 = theta*theta;
  theta3 = theta2*theta;
  if ( abs( theta ) > 0.000001 )
    a = (   1   - cos( theta ) ) / ( theta2 );
    b = ( theta - sin( theta ) ) / ( theta3 );
  else
    a = ( 1/2 ) * ( 0.5 - ( theta2/12 )*( 1 - ( theta2/30 )*( 1 - theta2/56 ) ) );
    b = ( 1/6 ) * (  1  - ( theta2/20 )*( 1 - ( theta2/42 )*( 1 - theta2/72 ) ) );
  end
  V = eye( 3 ) + a*wx + b*wx*wx;
  incr = [ [ expm( wx ),  V*u ]; [ 0 0 0 1 ] ];

  T  = incr * T;

end

printf( '----------\n' )
Txs = T*Xs;
txs = Txs( 1:2,: ) ./ Txs( 3,: );
e   = ( txs - ys )( : );
norms(iter+1) = norm(e);
T
plot( norms )
