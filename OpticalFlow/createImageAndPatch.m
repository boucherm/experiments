% Creates an image hopefully optical-flow friendly
function I = createImageAndPatch( kind, width, height )

  I_width  = width;
  I_height = height;

  %--------------------------------------------------------------------------------------------
  %--- Quadratic growth
  %--------------------------------------------------------------------------------------------
  if strcmp( kind, 'quadratic' ) || strcmp( kind, 'shell' )

    for ii=0:I_width-1
     Ix_quad( :, ii+1 ) = ii*ii*128/(I_width*I_width);
    end
    for ii=0:I_height-1
     Iy_quad( ii+1, : ) = ii*ii*128/(I_height*I_height);
    end

    I_quad = Ix_quad + Iy_quad;
    % ^ problematic when image gradient is evaluated where the slope is steep

  end

  %--------------------------------------------------------------------------------------------
  %--- Linear growth, axial
  %--------------------------------------------------------------------------------------------
  % Needs to evaluate the image gradient on the diagonal since outside of it the pseudo-hessian
  % is rank 1 only
  if strcmp( kind, 'linear' ) || strcmp( kind, 'shell' )

    for ii=0:I_width-1
      Ix_lin( :, ii+1 ) = ii*255/(I_width);
    end
    for ii=0:I_height-1
      Iy_lin( ii+1, : ) = ii*255/(I_height);
    end

    I_lin = ( Ix_lin + Iy_lin );

  end

  %--------------------------------------------------------------------------------------------
  %--- Linear growth, radial
  %--------------------------------------------------------------------------------------------
  % Works well for translation, probably not so much if a rotation is involved
  if strcmp( kind, 'radial' ) || strcmp( kind, 'shell' )

    for ii=0:I_width-1
    for jj=0:I_height-1
      I_rad( jj+1, ii+1 ) = 256*sqrt( ii*ii + jj*jj )/ sqrt( I_width*I_width + I_height*I_height );
    end
    end

  end

  %--------------------------------------------------------------------------------------------
  %--- Combination of axial and radial linear growths
  %--------------------------------------------------------------------------------------------
  %I = 0.5 * I_lin + 0.5 * I_rad ;
  %I = 0.5 * I_lin + 0.5 * ( 255*ones( I_height, I_width ) - I_rad );
  %I = 0.25 * I_lin + 0.5 * ( 255*ones( I_height, I_width ) - I_rad );
  if strcmp( kind, 'shell' )

    I_shell = 0.15 * I_lin + 0.5 * ( 255*ones( I_height, I_width ) - I_rad ); % looks like a sea-shell

  end

  %--------------------------------------------------------------------------------------------
  %--- (<1) exponent
  %--------------------------------------------------------------------------------------------
  % y = x^(1/n)
  % =>
  % n = ln(x) / ln(y)
  %
  % x ∈ [ a, b ]
  % y ∈ [ c, d ]
  % =>
  % n ∈ [ ln(a)/ln(d), ln(b)/ln(c) ] = [ n_min, n_max ]
  %
  % m ∈ [ 0, 255 ]
  % => m = f(n) = ( n - n_min ) / ( n_max - n_min ) * 255
  if strcmp( kind, 'smexp' )

    u_min = 1;
    u_max = I_width;
    v_min = 1;
    v_max = I_height;
    x_min = u_min;   % use xmin, xmax, ymin, ymax
    x_max = u_max;   % to "move" the curves
    y_min = v_min+1; % in the image
    y_max = v_max+1; %
    %x_min = 1;
    %x_max = 5;
    %y_min = 2;
    %y_max = 5;
    n_min = log( x_min ) / log( y_max );
    n_max = log( x_max ) / log( y_min );
    coeff = 255/( n_max - n_min );
    for u=1:I_width
    for v=1:I_height

      x = x_min + ( u - u_min ) / ( u_max - u_min ) * ( x_max - x_min );
      y = y_min + ( v - v_min ) / ( v_max - v_min ) * ( y_max - y_min );

      n = log( x ) / log( y );
      m = ( n - n_min ) * coeff;

      I_smexp( v, u ) = m;

    end
    end

    I = I_smexp;

  end

  %--------------------------------------------------------------------------------------------
  %--- Scales of ln()
  %--------------------------------------------------------------------------------------------
  % y = n.ln(x) <=> n = y/ln(x)
  %
  % x ∈ [ a, b ]
  % y ∈ [ c, d ]
  % =>
  % n ∈ [ c/ln(b), d/ln(a) ] = [ n_min, n_max ]
  %
  % m ∈ [ 0, 255 ]
  % => m = f(n) = ( n - n_min ) / ( n_max - n_min ) * 255
  if strcmp( kind, 'slog' )

    u_min = 1;
    u_max = I_width;
    v_min = 1;
    v_max = I_height+1;
    %x_min = u_min+1;
    %x_max = u_max+1;
    %y_min = v_min;
    %y_max = v_max;
    x_min = 1.1;
    x_max = 1.5;
    y_min = 1.1;
    y_max = 5;
    n_min = y_min / log( x_max );
    n_max = y_max / log( x_min );
    coeff = 255/( n_max - n_min );
    for u=1:I_width
    for v=1:I_height

      x = x_min + ( u - u_min ) / ( u_max - u_min ) * ( x_max - x_min );
      y = y_min + ( v - v_min ) / ( v_max - v_min ) * ( y_max - y_min );

      n = y / log( x );
      m = ( n - n_min ) * coeff;

      I_slog( v, u ) = m;
    end
    end

    I = I_slog;

  end

  %--------------------------------------------------------------------------------------------
  %--------------------------------------------------------------------------------------------
  %--------------------------------------------------------------------------------------------
  if strcmp( kind, 'quadratic' )
    I = I_quad;
  end
  if strcmp( kind, 'linear' )
    I = I_lin;
  end
  if strcmp( kind, 'radial' )
    I = I_rad;
  end
  if strcmp( kind, 'smexp' )
    I = I_smexp;
  end
  if strcmp( kind, 'slog' )
    I = I_slog;
  end
  if strcmp( kind, 'shell' )
    I = I_shell;
  end


end
