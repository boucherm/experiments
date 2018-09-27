% Select, and interpolate, an image from a set of coordinates
% IN
%   I  : image
%   xs : matrix, size 2xN, for each pixel of Iv associates coordinates to pick in I
%   Mv : binary mask, size Nx1, allows to ignore some pixels of Iv
% OUT
%   Iv : result image, size Nx1
%   Mv : updated mask, defining the validity/binary weight of each pixel, size Nx1
function [ Iv, Mv ] = mapImage( I, xs, Mv )

  N  = size( xs, 2 );
  Iv = zeros( N, 1 );

  h = size( I, 1 );
  w = size( I, 2 );

  Xs = xs( 1, : )';
  Ys = xs( 2, : )';

  Mv = Mv .* ( Xs >= 1 )       ...
          .* ( Xs <= w )       ...
          .* ( Ys >= 1 )       ...
          .* ( Ys <= h )       ...
          .* ( 1 - isnan(Xs) ) ...
          .* ( 1 - isnan(Ys) );

  Xsmin = floor( Xs );
  Ysmin = floor( Ys );
  Xsmax = ceil( Xs );
  Ysmax = ceil( Ys );

  Wsxmin = Xsmax - Xs;
  Wsymin = Ysmax - Ys;
  Wsxmax = 1 - Wsxmin;
  Wsymax = 1 - Wsymin;

  M = logical( Mv );

  idYsminXsmin = sub2ind( size(I), Ysmin(M), Xsmin(M) );
  idYsmaxXsmin = sub2ind( size(I), Ysmax(M), Xsmin(M) );
  idYsminXsmax = sub2ind( size(I), Ysmin(M), Xsmax(M) );
  idYsmaxXsmax = sub2ind( size(I), Ysmax(M), Xsmax(M) );

  Iv( M ) = Wsxmin(M) .* Wsymin(M) .* I( idYsminXsmin ) ...
          + Wsxmin(M) .* Wsymax(M) .* I( idYsmaxXsmin ) ...
          + Wsxmax(M) .* Wsymin(M) .* I( idYsminXsmax ) ...
          + Wsxmax(M) .* Wsymax(M) .* I( idYsmaxXsmax );

end


%% Unvectorized, slow, but understandable
%function [ Iv, Mv ] = mapImage( I, xs, Mv )
%
%  N  = size( xs, 2 );
%  Iv = zeros( N, 1 );
%
%  h = size( I, 1 );
%  w = size( I, 2 );
%
%  for ii=1:length( Iv )
%    if ( 0 == Mv( ii ) )
%      continue;
%    end
%
%    x = xs( 1, ii );
%    y = xs( 2, ii );
%
%    if ( ( x < 1 ) || ( y < 1 ) || ( x > w ) || ( y > h ) || isnan( x ) || isnan( y ) )
%      Mv( ii ) = 0;
%      continue;
%    else
%      Mv( ii ) = 1;
%    end
%
%    xmax =  ceil( x );
%    xmin = floor( x );
%
%    ymax =  ceil( y );
%    ymin = floor( y );
%
%    wxmin = xmax - x;
%    wxmax = 1 - wxmin;
%
%    wymin = ymax - y;
%    wymax = 1 - wymin;
%
%    %[ x, y, ymin, xmin, ymax, xmax ]
%
%    Iv( ii ) = wxmin * wymin * I( ymin, xmin ) ...
%             + wxmin * wymax * I( ymax, xmin ) ...
%             + wxmax * wymin * I( ymin, xmax ) ...
%             + wxmax * wymax * I( ymax, xmax );
%  end
%
%end
