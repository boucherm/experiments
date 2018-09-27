% Defines ( among other things ):
%  rI       : reference image
%  cI       : current   _____
%  rIv      : reference image layed out as a column vector ( as [ col1; col2; ...; colN ] )
%  fu       : camera intrinsics
%  fv       : ______ __________
%  K        : ______ __________ matrix
%  c_width  : camera image width  ( i.e. width  of rI and cI )
%  c_height : ______ _____ height _ ____ height __ __ ___ __ _
%  ps       : vector of coordinates of the pixels in an image
%  cTr      : transform from R to C
%  rD       : depth of image R
%  cD       : _____ __ _____ C


%---
% Generate texture
%I = createImageAndPatch( 'smexp', 256, 256 );
%I = createImageAndPatch( 'quadratic', 256, 256 );
%I = createImageAndPatch( 'shell', 256, 256 );
%I = createImageAndPatch( 'radial', 256, 256 );
%I = createImageAndPatch( 'slog', 256, 256 );
%I = double( imread( 'check_mono.png' ) );
I = double( imread( 'pattern.png' ) );
[ i_height, i_width ] = size( I );

%---
% Define a plane / 4 points
P1 = [ -1; -1;  0 ];
P2 = [ -1;  1;  0 ];
P3 = [  1;  1;  0 ];
P4 = [  1; -1;  0 ];
n  = [  0;  0; -1 ];

%---
% Define an image to plane similarity transform
pSi = [ 2/i_width , 0          , -1;
          0       , 2/i_height , -1;
          0       , 0          ,  1 ];
iSp = inv( pSi );

%---
% Define first pose
rRp = expm( hat( [ 0.1; 0.1; 0.0 ] ) );
rtp = [ 0.0; 0.2; 2.2 ];
rTp = [ [ rRp, rtp ]; [ 0 0 0 1 ] ];

%---
% Define second pose and relative pose
cRr = expm( hat( [ 0.1; 0.1; 0.2 ] ) );
ctr = [ -0.1; 0.3; -0.2 ];
cTr = [ [ cRr, ctr ]; [ 0 0 0 1 ] ];

cTp = cTr * rTp;

%---
% Define intrinsics
scale    = 2^4;
c_width  = 32 * scale;
c_height = 24 * scale;
fu       = 50 * scale;
fv       = 50 * scale;
cu       = c_width/2;
cv       = c_height/2;
K = [ fu , 0  , cu;
       0 , fv , cv;
       0 , 0  , 1 ];
Kinv = inv( K );

%---
% Compute plane-image homographies for each pose

% Use the fact that the plane is Z=0
% see Zhang's calibration section 2.2 for explanation
rHp = K*[ rTp( 1:3, [ 1 2 4 ] ) ];
cHp = K*[ cTp( 1:3, [ 1 2 4 ] ) ];
pHr = inv( rHp );
pHc = inv( cHp );
iHr = iSp * pHr;
iHc = iSp * pHc;

%---
% Generate images from homographies and texture

% Generate the vector of coordinates needed to scan the camera image
[ temp1, temp2 ] = meshgrid( 1:c_width, 1:c_height );
ps   = [ temp1(:)'; temp2(:)' ];
N    = c_width*c_height;
mask = ones( c_width*c_height, 1 );

% Reference
rIv = mapImage( I, warpSl3(ps,iHr), mask );
rI  = reshape( rIv, [c_height,c_width] );
% Current
cIv = mapImage( I, warpSl3(ps,iHc), mask );
cI  = reshape( cIv, [c_height,c_width] );

%---
% Generate depth maps for each image

% Reference
ys  = warpSl3( ps, pHr );
pYs = [ ys; zeros(1,N); ones(1,N) ];
rYs = rTp * pYs;
rD  = reshape( rYs( 3, : ), [c_height,c_width] );

% Current
ys  = warpSl3( ps, pHc );
pYs = [ ys; zeros(1,N); ones(1,N) ];
cYs = cTp * pYs;
cD  = reshape( cYs(3,:), [c_height,c_width] );


% Punch holes in rD
%floor(0*c_height/5)
%floor(1*c_height/5)
%floor(0*c_width/5)
%floor(1*c_width/5)
rD( floor(0*c_height/5)+1:floor(1*c_height/5), floor(0*c_width/5)+1:floor(1*c_width/5) ) = 0;
rD( floor(1*c_height/5)+1:floor(2*c_height/5), floor(4*c_width/5)+1:floor(5*c_width/5) ) = 0;
rD( floor(2*c_height/5)+1:floor(3*c_height/5), floor(2*c_width/5)+1:floor(3*c_width/5) ) = 0;
rD( floor(4*c_height/5)+1:floor(5*c_height/5), floor(1*c_width/5)+1:floor(2*c_width/5) ) = 0;

for i=0:2
for j=0:2
  rD( 4+i:10:end, 4+j:10:end ) = 0;
end
end


%---
% Drawings

%figure
%error_img = diffImagePatch( cI, rI );
%surf( 1:size(cI,2), 1:size(cI,1), error_img );

draw_generated = 0;
if ( draw_generated )
  figure;
  subplot( 3, 2, 1 );
  axis on;
  grid on;
  hold on;
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  Poses( :, :, 1 ) = rTp;
  Poses( :, :, 2 ) = cTp;
  drawCams( Poses );
  line1 = [ P1, P2 ];
  line2 = [ P2, P3 ];
  line3 = [ P3, P4 ];
  line4 = [ P4, P1 ];
  plot3( line1(1,:), line1(2,:), line1(3,:), 'color', [ 1, 0, 0 ] );
  plot3( line2(1,:), line2(2,:), line2(3,:), 'color', [ 0, 1, 0 ] );
  plot3( line3(1,:), line3(2,:), line3(3,:), 'color', [ 0, 0, 1 ] );
  plot3( line4(1,:), line4(2,:), line4(3,:), 'color', [ 0, 0, 0 ] );
  ;

  colormap( 'gray' )
  subplot( 3, 2, 2 ); imshow(  I/255 );

  subplot( 3, 2, 3 ); imshow( rI/255 );
  min_rd = min( min( rD ) );
  max_rd = max( max( rD ) );
  min_cd = min( min( cD ) );
  max_cd = max( max( cD ) );
  min_d  = min( [ min_rd, min_cd ] );
  max_d  = max( [ max_rd, max_cd ] );
  subplot( 3, 2, 4 ); imshow( ( rD - min_d ) / ( max_d - min_d ) );

  subplot( 3, 2, 5 ); imshow( cI/255 );
  subplot( 3, 2, 6 ); imshow( ( cD - min_d ) / ( max_d - min_d ) );
end
