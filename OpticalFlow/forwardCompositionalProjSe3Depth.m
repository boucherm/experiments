% Following "Lucas-Kanade: 20 years on"
%   . reference plays the "template"
%   . current   plays the "image"
%   . the warp maps coordinates from template to image

function cTr = forwardCompositionalSe3( rI, cI, rD, K, cTr, n_iter, step, draw )
  pkg load image;

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
  % (4)
  %tic;
  xs = [ ps; ones(1,N) ];
  xs = Kinv * xs;
  Xs = xs .* ds;
  Xs = [ Xs; ones(1,N) ];

  [ G_yz, G_zx, G_xy, G_x, G_y, G_z ] = se3Generators();

  Gyz_Xs = G_yz * Xs; % R_yz
  Gzx_Xs = G_zx * Xs; % R_zx
  Gxy_Xs = G_xy * Xs; % R_xy
  Gx_Xs  = G_x  * Xs; % t_x
  Gy_Xs  = G_y  * Xs; % t_y
  Gz_Xs  = G_z  * Xs; % t_z
  %'se3 jacobian precomputation'
  %toc

  % Unvectorized, slow, but understandable
  %tic
  %for ii=1:N
  %  if 0 == Xs( 3, ii )
  %    continue
  %  end
  %  Jproj = 1/Xs( 3, ii ) * [ [ fu , 0  , -fu*Xs( 1 , ii )/Xs( 3 , ii ) , 0 ]; ...
  %                            [  0 , fv , -fv*Xs( 2 , ii )/Xs( 3 , ii ) , 0 ] ];

  %  i_one = ii*2 - 1;
  %  i_snd = ii*2;

  %  Jw( i_one:i_snd, 1 ) = Jproj * Gyz_Xs( :, ii );
  %  Jw( i_one:i_snd, 2 ) = Jproj * Gzx_Xs( :, ii );
  %  Jw( i_one:i_snd, 3 ) = Jproj * Gxy_Xs( :, ii );
  %  Jw( i_one:i_snd, 4 ) = Jproj *  Gx_Xs( :, ii );
  %  Jw( i_one:i_snd, 5 ) = Jproj *  Gy_Xs( :, ii );
  %  Jw( i_one:i_snd, 6 ) = Jproj *  Gz_Xs( :, ii );
  %end
  %'warp jacobian precomputation'
  %toc

  %tic
  Jproj = repmat( [ [ fu, 0, -fu, 0 ]; ...
                    [ 0, fv, -fv, 0 ] ], N, 1 );
  Xs12_v = reshape( Xs(1:2,:), 2*N, 1 );
  Xs33_v = reshape( [Xs(3,:);Xs(3,:)], 2*N, 1 );
  Jproj( :, 1 ) = Jproj(:,1) ./ Xs33_v;
  Jproj( :, 2 ) = Jproj(:,2) ./ Xs33_v;
  Jproj( :, 3 ) = Jproj(:,3) .* Xs12_v ./ ( Xs33_v .* Xs33_v );
  % Jproj( :, 4 ) is all zeros
  Jproj( isnan(Jproj) ) = 0;
  Jproj( isinf(Jproj) ) = 0;
  %'projection jacobian precomputation'
  %toc

  %tic
  Jw        = zeros( 2*N, 6 );
  u_indices = 1:2:2*N;
  v_indices = 2:2:2*N;
  Jw( u_indices, 1 ) = sum( Jproj(u_indices,:)' .* Gyz_Xs );
  Jw( v_indices, 1 ) = sum( Jproj(v_indices,:)' .* Gyz_Xs );
  Jw( u_indices, 2 ) = sum( Jproj(u_indices,:)' .* Gzx_Xs );
  Jw( v_indices, 2 ) = sum( Jproj(v_indices,:)' .* Gzx_Xs );
  Jw( u_indices, 3 ) = sum( Jproj(u_indices,:)' .* Gxy_Xs );
  Jw( v_indices, 3 ) = sum( Jproj(v_indices,:)' .* Gxy_Xs );
  Jw( u_indices, 4 ) = sum( Jproj(u_indices,:)' .*  Gx_Xs );
  Jw( v_indices, 4 ) = sum( Jproj(v_indices,:)' .*  Gx_Xs );
  Jw( u_indices, 5 ) = sum( Jproj(u_indices,:)' .*  Gy_Xs );
  Jw( v_indices, 5 ) = sum( Jproj(v_indices,:)' .*  Gy_Xs );
  Jw( u_indices, 6 ) = sum( Jproj(u_indices,:)' .*  Gz_Xs );
  Jw( v_indices, 6 ) = sum( Jproj(v_indices,:)' .*  Gz_Xs );
  %'warp jacobian precomputation'
  %toc

  %---------------------------------
  %fprintf( '###################\n' )
  errors = zeros( 1, n_iter+1 );
  for iter=1:n_iter+1

    %---------------------------------
    % (1)
    %tic;
    qs = warpProjSe3Depth( ps, ds, cTr, K );
    %'coordinates warping'
    %toc
    %tic;
    Mv          = ( ds > 0 )';
    [ Iwv, Mv ] = mapImage( cI, qs, Mv );
    Iw          = reshape( Iwv, size(rI) );
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
    %fprintf('--- %d\n', iter)
    %fprintf( 'error: %f\n', norm( e_cur ) );
    errors( iter ) = e_cur;

    if draw
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
    M       = reshape( Mv, size(rI) );
    Mw      = imerode( M, [ 1, 1, 1 ] ); % need to ignore gradient at the boundaries of the warped
    Mh      = imerode( M, [ 1; 1; 1 ] ); % image validity mask
    Jiw_x   = Mw .* Xdiff( Iw );
    Jiw_y   = Mh .* Ydiff( Iw );
    Jiw_x_v = reshape( Jiw_x, c_height*c_width, 1 );
    Jiw_y_v = reshape( Jiw_y, c_height*c_width, 1 );
    Jiw_v   = [ Jiw_x_v, Jiw_y_v ];
    %'reference image gradient computation'
    %toc

    %---------------------------------
    % (5)
    %tic;

    % Unvectorized, slow, but understandable
    %for ii=1:N
    %  i_one = ii*2 - 1;
    %  i_snd = ii*2;
    %  J( ii, : ) = Jiw_v( ii, : ) * Jw( i_one:i_snd, : );
    %end

    % Vectorized
    J      = zeros( N, 6 );
    Jiw_v2 = reshape( Jiw_v', 2*N, 1 );
    for e=1:6
      temp      = Jiw_v2 .* Jw( :, e );
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
    delta = -pinv(H) * Jt_e; % or Levenberg-Marquardt
    delta = step * delta;

    %---------------------------------
    % (9)
    incr = expSe3( delta );
    cTr  = cTr * incr;
    %'solve and update'
    %toc

  end
  fprintf( 'Performed %d iterations\n', iter );
end
