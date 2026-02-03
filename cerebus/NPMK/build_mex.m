function build_mex()
% BUILD_MEX Compiles the readNSxMex C++ function
% Detects platform and sets flags

srcFile = 'readNSxMex.cpp';
outputName = 'readNSxMex';

fprintf('Compiling %s...\n', srcFile);

if ispc
    % Windows
    % Often requires MinGW or VS.
    % We assume user has "mex -setup" configured.
    mex(srcFile, '-largeArrayDims', '-output', outputName, '-O', 'COMPFLAGS="$COMPFLAGS /O2"');
elseif isunix
    % Linux/Mac
    mex(srcFile, '-largeArrayDims', '-output', outputName, '-O', 'CXXFLAGS="$CXXFLAGS -O3"');
else
    mex(srcFile, '-largeArrayDims', '-output', outputName, '-O');
end

fprintf('Compilation successful. Created %s.%s\n', outputName, mexext);

end
