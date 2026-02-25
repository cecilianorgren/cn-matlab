% ADS database
%
save_url ='/Users/Cecilia/cn-matlab/ads/';
save_filename = 'ads_get.xml';

url = 'http://adsabs.harvard.edu/cgi-bin/';
action = 'abs_connect?';
bibcode='bibcode=2005PhPl...12j2110C';
query_type='query_type=CITES';
data_type='data_type=PLAINTEXT';
data_type='data_type=XML';

[file,status]=urlwrite([url,action,bibcode,'&',query_type,'&',data_type],[save_url,save_filename]);
%open([save_url,save_filename])
%%
xDoc = xmlread([save_url,save_filename]);
xRoot=xDoc.getDocumentElement;
hej=char(xRoot.getAttribute('bibcode'));
%%
z = xmltools([save_url,save_filename]);