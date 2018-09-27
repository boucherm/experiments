function [ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators()

  G_yz = [ [ hat( [ 1; 0; 0 ] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];
  G_zx = [ [ hat( [ 0; 1; 0 ] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];
  G_xy = [ [ hat( [ 0; 0; 1 ] ), zeros( 3, 1 ) ]; [ zeros( 1, 4 ) ] ];
  G_x  = zeros( 4 ); G_x( 1, 4 ) = 1;
  G_y  = zeros( 4 ); G_y( 2, 4 ) = 1;
  G_z  = zeros( 4 ); G_z( 3, 4 ) = 1;

end
