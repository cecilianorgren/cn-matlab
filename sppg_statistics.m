imember = 0;
imember = imember + 1;
group(imember).author = 'Norgren, Cecilia';
group(imember).year = '2017-2022';
imember = imember + 1;
group(imember).author = 'Tenfjord, Paul';
group(imember).year = '2017-2020';
imember = imember + 1;
group(imember).author = 'Hesse, Michael';
group(imember).year = '2016-2019';



search_string = '';
for im = 1:numel(group)
  search_string = [search_string,sprintf('(author:"%s" AND year:%s)',group(im).author,group(im).year),' OR '];
end

search_string = search_string(1:end-4);
disp(search_string)
