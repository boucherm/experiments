% Following "Lucas-Kanade: 20 years on"
%   . reference plays the "template"
%   . current   plays the "image"
%   . the warp maps coordinates from template to image

function cTr = forwardAdditiveSe3( rI, cI, rD, K, cTr, n_iter, step, draw )

  %---------------------------------
  if 8 > nargin
    draw = false;
  end

  if draw
    figure;
    subplot( 3, 2, 1 ); imshow( rI/255 );
    title( 'ref image' );
    subplot( 3, 2, 2 ); imshow( cI/255 );
    title( 'cur image' );
  end

  %---------------------------------
  %tic;
  r_height = size( rI, 1 );
  r_width  = size( rI, 2 );
  c_height = size( cI, 1 );
  c_width  = size( cI, 2 );

  Kinv = inv( K );
  fu   = K( 1, 1 );
  fv   = K( 2, 2 );

  rIv = double( reshape( rI, r_height*r_width, 1 ) );

  [ temp1, temp2 ] = meshgrid( 1:r_width, 1:r_height );
  ps               = [ temp1(:)'; temp2(:)' ];
  ds               = rD( sub2ind( size(rD), ps(2,:), ps(1,:) ) );
  N                = size( ps, 2 );

  %'preparation'
  %toc

  %---------------------------------
  % (3) - prepare
  %tic;
  Ji_x = Xdiff( cI );
  Ji_y = Ydiff( cI );
  %'current image gradient computation'
  %toc

  %---------------------------------
  % (4) - prepare
  %tic;
  xs = [ ps; ones(1,N) ];
  xs = Kinv * xs;
  Xs = xs .* ds;
  Xs = [ Xs; ones(1,N) ];

  [ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();
  %'se3 jacobian preparation'
  %toc

  %---------------------------------
  %fprintf( '###################\n' )
  errors = zeros( 1, n_iter+1 );
  for iter=1:n_iter+1

    %---------------------------------
    % (1)
    %tic;
    qs          = warpProjSe3Depth( ps, ds, cTr, K );
    %'coordinates warping'
    %toc
    Mv          = ( ds > 0 )';
    [ Iwv, Mv ] = mapImage( cI, qs, Mv );
    %'image mapping'
    %toc

    %---------------------------------
    % (2)
    %tic;
    residual = ( rIv - Iwv );
    residual = residual .* Mv;
    e_cur    = residual' * residual;
    %'residuals evaluation'
    %toc

    %---------------------------------
    % error control and drawings
    %fprintf( '--- %d\n', iter )
    %fprintf( 'error: %f\n', norm( e ) );
    errors( iter ) = e_cur;

    if draw
      Iw = reshape( Iwv, size( rI ) );
      subplot( 3, 2, 3 ); imshow( Iw/255 );
      title( 'cur warped as ref' );
      D = abs( Iw - rI );
      if 1 == iter
        m = max( max( D ) );
        if m > 0
          d_reg = 1 / m;
        end
      end
      subplot( 3, 2, 4 ); imshow( D*d_reg .* reshape( Mv, size(Iw) ) );
      title( 'diff' );
      subplot( 3, 2, [5,6] ); semilogy( errors(1:iter) );
      title( 'errors' );
      drawnow;
    end

    if ( n_iter+1 == iter ) || ( iter > 1 && e_cur > 0.99*e_prev )
      if ( e_cur > e_prev )
        'error increase'
        cTr = cTr_prev;
      end
      break;
    end

    e_prev   = e_cur;
    cTr_prev = cTr;

    %---------------------------------
    % (3)
    %tic;
    Jix_w_v = mapImage( Ji_x, qs, Mv );
    Jiy_w_v = mapImage( Ji_y, qs, Mv );
    Ji_w_v  = [ Jix_w_v, Jiy_w_v ];
    %'current image gradient warping'
    %toc

    %---------------------------------
    % (4)
    %tic;
    Ys     = cTr * Xs;
    Gyz_Ys = G_yz * Ys; % R_yz
    Gzx_Ys = G_zx * Ys; % R_zx
    Gxy_Ys = G_xy * Ys; % R_xy
    Gx_Ys  = G_x  * Ys; % t_x
    Gy_Ys  = G_y  * Ys; % t_y
    Gz_Ys  = G_z  * Ys; % t_z

    % Unvectorized, slow, but understandable
    %for ii=1:N
    %  Jproj = 1/Ys( 3, ii ) * [ [ fu , 0  , -fu*Ys( 1 , ii )/Ys( 3 , ii ) , 0 ]; ...
    %                            [  0 , fv , -fv*Ys( 2 , ii )/Ys( 3 , ii ) , 0 ] ];

    %  i_one = ii*2 - 1;
    %  i_snd = ii*2;

    %  Jw( i_one:i_snd, 1 ) = Jproj * Gyz_Ys( :, ii );
    %  Jw( i_one:i_snd, 2 ) = Jproj * Gzx_Ys( :, ii );
    %  Jw( i_one:i_snd, 3 ) = Jproj * Gxy_Ys( :, ii );
    %  Jw( i_one:i_snd, 4 ) = Jproj *  Gx_Ys( :, ii );
    %  Jw( i_one:i_snd, 5 ) = Jproj *  Gy_Ys( :, ii );
    %  Jw( i_one:i_snd, 6 ) = Jproj *  Gz_Ys( :, ii );
    %end

    % Vectorized
    Jproj = repmat( [ [ fu, 0, -fu, 0 ]; ...
                      [ 0, fv, -fv, 0 ] ], N, 1 );
    Ys12_v = reshape( Ys(1:2,:), 2*N, 1 );
    Ys33_v = reshape( [Ys(3,:); Ys(3,:)], 2*N, 1 );
    Jproj( :, 1 ) = Jproj( :, 1 ) ./ Ys33_v;
    Jproj( :, 2 ) = Jproj( :, 2 ) ./ Ys33_v;
    Jproj( :, 3 ) = Jproj( :, 3 ) .* Ys12_v ./ ( Ys33_v .* Ys33_v );
    % Jproj( :, 4 ) is all zeros
    Jproj( isnan( Jproj ) ) = 0;
    Jproj( isinf( Jproj ) ) = 0;
    %'projection jacobian computation'
    %toc

    %tic
    u_indices = 1:2:2*N;
    v_indices = 2:2:2*N;
    Jw        = zeros( 2*N, 6 );
    Jw( u_indices, 1 ) = sum( Jproj(u_indices,:)' .* Gyz_Ys );
    Jw( v_indices, 1 ) = sum( Jproj(v_indices,:)' .* Gyz_Ys );
    Jw( u_indices, 2 ) = sum( Jproj(u_indices,:)' .* Gzx_Ys );
    Jw( v_indices, 2 ) = sum( Jproj(v_indices,:)' .* Gzx_Ys );
    Jw( u_indices, 3 ) = sum( Jproj(u_indices,:)' .* Gxy_Ys );
    Jw( v_indices, 3 ) = sum( Jproj(v_indices,:)' .* Gxy_Ys );
    Jw( u_indices, 4 ) = sum( Jproj(u_indices,:)' .*  Gx_Ys );
    Jw( v_indices, 4 ) = sum( Jproj(v_indices,:)' .*  Gx_Ys );
    Jw( u_indices, 5 ) = sum( Jproj(u_indices,:)' .*  Gy_Ys );
    Jw( v_indices, 5 ) = sum( Jproj(v_indices,:)' .*  Gy_Ys );
    Jw( u_indices, 6 ) = sum( Jproj(u_indices,:)' .*  Gz_Ys );
    Jw( v_indices, 6 ) = sum( Jproj(v_indices,:)' .*  Gz_Ys );
    %'warp jacobian computation'
    %toc


    %Jw
    %---------------------------------
    % (5)
    %tic;

    % Unvectorized, slow, but understandable
    %for ii=1:N
    %  i_one = ii*2 - 1;
    %  i_snd = ii*2;

    %  J( ii, : ) = Ji_w_v( ii, : ) * Jw( i_one:i_snd, : );
    %end

    % Vectorized
    J       = zeros( N, 6 );
    Ji_w_v2 = reshape( Ji_w_v', 2*N, 1 );
    for e=1:6
      temp      = Ji_w_v2 .* Jw( :, e );
      J( :, e ) = temp( 1:2:end ) + temp( 2:2:end );
    end

    %'jacobian computation'
    %toc

    %---------------------------------
    % (6)
    %tic;
    H = J' * J;
    %'hessian computation'
    %toc

    %---------------------------------
    % (7)
    %tic;
    Jt_e = J' * residual;
    %'Jte'
    %toc

    %tic;
    %---------------------------------
    % (8)
    delta = -pinv(H) * Jt_e;
    delta = step * delta;

    %---------------------------------
    % (9)
    incr = expSe3( delta );
    cTr  = incr * cTr;
    %'solve and update'
    %toc

  end
  fprintf( 'Performed %d iterations\n', iter );
end
