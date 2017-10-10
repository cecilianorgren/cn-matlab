% finn.no

%resultssearch = urlread('https://www.finn.no/realestate/homes/search.html?published=1&q=bergen');
%resultssearch = urlread('https://www.finn.no/realestate/homes/search.html?filters=');
% str_base = 'https://www.finn.no/realestate/homes/search.html?';
% str_published = 'published=1';
% str_location = 'q=bergen';
% str_search = [str_base str_published str_location];
% results_search = urlread(str_search);
latitude_office = 60.3834109;
longitude_office = 5.3286653;

str_base = 'https://www.finn.no/realestate/homes/search.json?';
str_published = '';'published=';
str_location = 'q=bergen';

str_search = [str_base str_published str_location];

json_results = webread(str_search);
objects = json_results.displaySearchResults;
%%
units = irf_units;
for ii = 1:numel(objects)
  if isfield(objects{ii},'bodyRow')
    result_object = webread(['https://www.finn.no' objects{ii}.adUrl]);
    cell_coord = regexp(result_object,'&lat=(\d*.\d*)|&lng=(\d*.\d*)','tokens');
    latitude = str2num(cell_coord{1}{1});
    longitude = str2num(cell_coord{2}{1});
    distance = acos( sind(latitude_office)*sind(latitude) + cosd(latitude_office)*cosd(latitude)*cosd(longitude_office-longitude))*units.RE*1e-3;
    size = regexp(objects{ii}.bodyRow{1},'(\d*)','match');
    price = regexp(objects{ii}.bodyRow{2},'\s-\s|,-','split'); price = price(1:end-1);    
    price_ = price{1}; price_(isspace(price_))=[];
    price_m2 = str2num(price_)/str2num(size{1});
    %disp(sprintf('%3.0f: %7.4f, %7.4f, %3.1f km , %6s, %11s, %s',ii,latitude,longitude,distance,objects{ii}.bodyRow{1},objects{ii}.bodyRow{2},objects{ii}.topRowCenter))
    if distance<10 && str2num(size{1})<100;
      disp(sprintf('%3.0f: %3.1f km, %.0f kNOK/m2, %6s, %11s, %s',ii,distance,price_m2*1e-3,objects{ii}.bodyRow{1},objects{ii}.bodyRow{2},objects{ii}.topRowCenter))    
    else
      %disp(sprintf('%3.0f:',ii))
    end
  else
    %disp(sprintf('%3.0f:',ii))
  end
end
%% 
result_object = webread('https://www.finn.no/realestate/homes/ad.html?finnkode=93490622');

%%
results_object = urlread('https://www.finn.no/realestate/homes/ad.html?finnkode=93490622');
loc_address =  findstr(results_object,'Hulda Garborgs gate');

for ii = 1:numel(loc_address);
  results_object(loc_address(ii)+[-1000:1000]) 
end


%% Coordinates
cell_coord = regexp(results_object,'&lat=(\d*.\d*)|&lng=(\d*.\d*)','tokens');
latitude = str2num(cell_coord{1}{1});
longitude = str2num(cell_coord{2}{1});

%% Address
cell_coord = regexp(results_object,'&lat=(\d*.\d*)|&lng=(\d*.\d*)','tokens');
latitude = str2num(cell_coord{1}{1});
longitude = str2num(cell_coord{2}{1});

%% Get direction from Google
%https://maps.googleapis.com/maps/api/geocode/json?latlng={0}&key={1}
%'https://maps.googleapis.com/maps/api/directions/json?origin=75+9th+Ave+New+York,+NY&destination=MetLife+Stadium+1+MetLife+Stadium+Dr+East+Rutherford,+NJ+07073&key=YOUR_API_KEY'
