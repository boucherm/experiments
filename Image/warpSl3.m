% SL(3) warp of coordinates
% IN
%   ps : coordinates to warp ( should be size 2xN or 3xN for N coordinates )
%    T : SL(3) matrix transformation
% OUT
%   qs : warped coordinates

function qs = warpSl3( ps, H )

  N = size( ps, 2 );

  % make sure ps is size 3 x N
  if ( 2 == size( ps, 1 ) )
    ps = [ ps; ones( 1, N ) ];
  end

  qs = H * ps;

  qs = qs( 1:2, : ) ./ qs( 3, : );

end
