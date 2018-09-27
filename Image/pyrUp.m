% TODO DOC
function pyr = pyrUp( I, levels, M )

  pkg load image;

  masked = false;
  if nargin > 2
    masked = true;
  end

  [ h, w ] = size( I );
  pyr      = cell( levels, 1 );
  pyr{1}   = I;

  for lvl = 2:levels
    if ! masked
      pyr{lvl} = imReduce( pyr{lvl-1} );
    else

      h     = h/2;
      w     = w/2;
      Iprev = pyr{ lvl-1 };
      Mprev = M;
      I     = zeros( h, w );
      M     = zeros( h, w );

      % Unvectorized, slow, but understandable
      %for u = 1:w
      %for v = 1:h
      %  v0 = 2*v-1; v1 = 2*v;
      %  u0 = 2*u-1; u1 = 2*u;
      %  s  = Mprev(v0,u0) + Mprev(v1,u0) + Mprev(v0,u1) + Mprev(v1,u1);
      %  if s > 0
      %    I( v, u ) = 1/s * ( Iprev(v0,u0) * Mprev(v0,u0)
      %                      + Iprev(v1,u0) * Mprev(v1,u0)
      %                      + Iprev(v0,u1) * Mprev(v0,u1)
      %                      + Iprev(v1,u1) * Mprev(v1,u1) );
      %  end
      %  M( v, u ) = max( max( Mprev( v0:v1, u0:u1 ) ) );
      %end
      %end
      %pyr{lvl} = I;

      % Evaluate regularizer of each destination pixel
      M00 = Mprev;
      M01 = [ M00(2:end,   :  ); M00(end,  : ) ];
      M10 = [ M00(  :  , 2:end), M00( : , end) ];
      M11 = [ M10(2:end,   :  ); M10(end,  : ) ];

      S     = ( M00 + M01 + M10 + M11 )( 1:2:end, 1:2:end );
      Sm    = S > 0;
      R     = zeros( h, w );
      R(Sm) = 1./S(Sm);
      % To understand the next instruction, try this one:
      % A=magic(4), k=ones(2), C=conv2(A,k), S=C(2:2:end,2:2:end)
      I        = R .* ( conv2( Iprev .* Mprev, ones(2) )( 2:2:end, 2:2:end ) );
      pyr{lvl} = I;

      % Update mask
      Ms        = zeros( size(Mprev,1), size(Mprev,2), 4 );
      Ms(:,:,1) = M00;
      Ms(:,:,2) = M01;
      Ms(:,:,3) = M10;
      Ms(:,:,4) = M11;
      M         = max( Ms, [], 3 )( 1:2:end, 1:2:end );

    end
  end

end


function J = imReduce( I )

  Gx = 0.25 * [ 1, 2, 1 ];
  Gy = 0.25 * [ 1; 2; 1 ];

  I = conv2( [ I(:,1), I, I( : ,end) ], Gx )(   :   ,3:end-2);
  I = conv2( [ I(1,:); I; I(end, : ) ], Gy )(3:end-2,   :   );
  J = I( 1:2:end, 1:2:end );

end
