function listgpus()
% listgpus()
%
% List out as much information about the GPUs in the computer as possible.
% current it only dumps the content of the cudaGetDeviceProperties call.

functions = {'listgpus_cuda50', 'listgpus_cuda42', 'listgpus_cuda40'};

for ii = functions
    try
        feval(ii{1});
        return;
    catch ME
        % don't do anything here, just supress the fails
    end
end

error(['Sorry, all versions of listgpus failed to run on your system. ' ...
    'Please verify Cuda is installed properly']);
end