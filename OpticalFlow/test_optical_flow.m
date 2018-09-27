more off;
clear all;
close all;

generate_data_proj_se3_depth;

cTr_   = eye(4);
n_iter = 20;
step   = 1.2;
%levels = 6;
levels = 5;
%levels = 4;
%levels = 3;
%levels = 2;

draw   = true;
%draw   = false;

rI_pyr = pyrUp( rI, levels );
cI_pyr = pyrUp( cI, levels );
rD_pyr = pyrUp( rD, levels, rD > 0 );

id = tic;
for lvl = levels:-1:1
  scale = 1.0 / power( 2, lvl-1 );
  Klvl  = [ K( 1:2, 1:3 ) * scale; 0 0 1 ];
  fprintf( '----------------- level: %d\n', lvl );
  fprintf( 'Intrinsics:\n' )
  Klvl
  fprintf( 'Input pose:\n' )
  cTr_
  %cTr_ = forwardAdditiveProjSe3Depth(
  %cTr_ = forwardCompositionalProjSe3Depth(
  cTr_ = inverseCompositionalProjSe3Depth(
           rI_pyr{lvl},
           cI_pyr{lvl},
           rD_pyr{lvl},
           Klvl,
           cTr_,
           n_iter,
           step,
           draw );
  fprintf( 'Result pose:\n' ); cTr_
  fprintf( 'Real pose:\n' );  cTr
  fprintf( 'Errors:\n' );
  E = inv( cTr ) * cTr_;
  [ vee( logm( E(1:3,1:3) ) )', E(1:3,4)' ]
  %keyboard
  %break
end
toc( id )
