function data1 = cn_offset(sc1,data1,cosys,offset)

% CN_OFFSET Calculate offset along spin axis for FGM-data.
% 
% Example: 
%   [offset, newB3] = cn_offset(sc3,B3,sc4,B4,cosys);
%   
%   where B3 B4 are time series, sc3 sc4 are the numbers 3 and 4
%   and  cosys is the coordinate system of B3 and B4, eg. cosys='GSM'.
%

if strcmp(lower(cosys),'dsi')
    data1(:,4)=data1(:,4)-offset;
else
    didata1=c_coord_trans(cosys,'dsi',data1,'cl_id',sc1);
    didata1(:,4)=didata1(:,4)-offset;
    data1=c_coord_trans('dsi',cosys,didata1,'cl_id',sc1);
end
