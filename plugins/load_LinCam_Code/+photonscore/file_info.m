function info = file_info(filename)
r = photonscore_mex(photonscore.Const.FN_FILE_INFO, filename);

info.photons_count = nan;
info.raw_info_ = r;

for i=1:size(r, 1)
  switch r{i, 1}
  case '/photons/TimerFrequency'
    info.aat_frequencty = str2num(r{i, 3});
  case 'Comment'
    info.comment = r{i, 3};
  case 'Created'
    info.created = r{i, 3};
  case '/photons/TacBias'
    info.dt_bias = str2num(r{i, 3});
  case '/photons/TacChannel'
    info.dt_channel = str2num(r{i, 3});
  case '/photons/DetectorGuid'
    info.detector_guid = r{i, 3};
  case 'FileGuid'
    info.file_guid = r{i, 3};
  case {'/photons/x', '/photons/y', '/photons/dt', '/photons/aat'}
    if r{i, 3} > 0
      info.photons_count = min(info.photons_count, r{i, 3});
    end
  case '/photons/ms'
    info.duration = double(max(r{i, 3})) / 1000;
  end
end

