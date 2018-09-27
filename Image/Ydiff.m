% Computes dI/dy
% Mirrors image at boundaries
function Ji_y = Ydiff( I )

  mask        = [ -1; 0; 1 ];
  Iextended   = [ I( 2, : ); I; I( end-1, : ) ];
  Jextended_y = conv2( Iextended, mask );
  Ji_y        = Jextended_y( 3:end-2, : );

end
