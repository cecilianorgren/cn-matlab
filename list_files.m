FL = struct([]);
for iYear = 2017
  for iMonth = 1:12
    FL_tmp = dir(sprintf('/Volumes/mms/mms3/fgm/brst/l2/%g/%02.0f/06/*.cdf',2017,6));
    FL = cat(1,FL,FL_tmp);
  end
end