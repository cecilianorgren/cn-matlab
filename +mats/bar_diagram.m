x = [2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015];
y = [12 4 14 28 34 29 28 36 27 31 46 42 43 42 49];


hca = subplot(1,1,1);
hb = bar(hca,x,y);
hb.FaceColor = [0 0.4470 0.7410];
hca.XTick = hca.XTick(1:2:end);
%hca.YLabel.String = 'Number of publications';
hca.Title.String = 'Number of publications';

