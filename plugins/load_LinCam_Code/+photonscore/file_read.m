function data = file_read(filename, dataset, offset, count)
if nargin ~= 4
  offset = 0;
  count = 0;

  info = photonscore.file_info(filename);
  for i=1:size(info.raw_info_, 1)
    if strcmp(info.raw_info_{i, 1}, dataset)
      count = info.raw_info_{i, 3}(1);
    end
  end

  if count == 0
    error 'Cannot find dataset';
  end
end

data = photonscore_mex(photonscore.Const.FN_FILE_READ, ...
    filename, dataset, offset, count);
