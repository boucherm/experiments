% Computes dI/dx
% Mirrors image at boundaries
function Ji_x = Xdiff( I )

  mask        = [ -1 0 1 ];
  Iextended   = [ I( :, 2 ), I, I( :, end-1 ) ];
  Jextended_x = conv2( Iextended, mask );
  Ji_x        = Jextended_x( :, 3:end-2 );

end
