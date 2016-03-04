% Script for compilation of TNT MEX files.
% Files are automatically shifted to the right directories for TNT to work.

try %% NearestNeighbor tracker
    if ispc % Windows with Microsoft Visual Studio compiler
        %% NearestNeighbor tracker
        run('NearestNeighborTracker\build_nn_tracker');
        copyfile('NearestNeighborTracker\mx_nn_tracker.mexw64','..\subfun\');
    elseif isunix % Unix (assuming gcc > 4.5)
        run('NearestNeighborTracker/build_nn_tracker');
        copyfile('NearestNeighborTracker/mx_nn_tracker.mexa64','../subfun/');
    end
catch err
    warning('Compilation of NearestNeighborTracker failed. ERROR: %s',err.message);
end

try %% TNT fitter
    if ispc % Windows with Microsoft Visual Studio compiler
        run('FastPsfFitting\build_psfFit_Image');
        copyfile('FastPsfFitting\ceres.dll','..\subfun\');
        copyfile('FastPsfFitting\mx_psfFit_Image.mexw64','..\subfun\');
    elseif isunix % Unix (assuming gcc > 4.5)
        run('FastPsfFitting/build_psfFit_Image');
        !cp TNTfitter/libceres* ../subfun/
        copyfile('FastPsfFitting/mx_psfFit_Image.mexa64','../subfun/');
    end
catch err
    warning('Compilation of FastPsfFitting (TNT fitter) failed. Did you build the ceres library first? -> see "COMPILATION_README.txt"\n ERROR: %s',err.message);
end