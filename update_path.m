% https://octave.org/doc/v4.0.1/Manipulating-the-Load-Path.html

% The following doesn't seem to work
%addpath( genpath( pwd, skip, [ pwd '/.git' ] ) );

addpath( [ pwd '/Image' ] );
addpath( [ pwd '/Lie' ] );
addpath( [ pwd '/OpticalFlow' ] );
addpath( [ pwd '/Optim' ] );
savepath();
