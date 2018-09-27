% SE(3)+projection warp of coordinates
% IN
%   ps : coordinates to warp ( should be size 2xN or 3xN for N coordinates )
%   ds : associated depths ( should be size 1xN ___ _ ___________ )
%    T : SE(3) matrix transformation
%    K : intrinsics parameters matrix
% OUT
%   qs : warped coordinates
function qs = warpProjSe3Depth( ps, ds, T, K )

  N    = size( ps, 2 );
  Kinv = inv( K );

  % make sure ps is size 3 x N
  if ( 2 == size( ps, 1 ) )
    ps = [ ps; ones( 1, N ) ];
  end

  xs = Kinv * ps;

  Xs = xs .* ds;
  Xs = [ Xs; ones( 1, N ) ];

  Ys = T * Xs;

  ys = Ys( 1:3, : ) ./ Ys( 3, : );
  qs = K * ys;
  qs = qs( 1:2, : );

end
