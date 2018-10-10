function [ G_yz, G_zx, G_xy ] = so3Generators()

  G_yz = hat( [ 1; 0; 0 ] );
  G_zx = hat( [ 0; 1; 0 ] );
  G_xy = hat( [ 0; 0; 1 ] );

end
