clear mex % so that we can overwrite the mexw64/mexa64 files

disp('Building NearestNeighborTracker ..')

if ispc % Windows with Microsoft Visual Studio compiler
    % Note: The /MT option statically links in the runtime, so it does not need to be deployed with the application.
    mex OPTIMFLAGS="/O2 /Oy- /DNDEBUG /MT" ...
        mx_nn_tracker.cpp   
elseif isunix % Unix (assuming gcc > 4.5)
    mex CXXFLAGS="\$CXXFLAGS -std=c++11 -O2 -static-libgcc -static-libstdc++" ...
        mx_nn_tracker.cpp
end