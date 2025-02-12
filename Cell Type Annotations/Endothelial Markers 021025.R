#https://www.biocompare.com/Editorial-Articles/598462-A-Guide-to-Endothelial-Cell-Markers/
DotPlot(obj, c('flt1',
               'cldn5a',
               'vwf',
               'ace',
               'cd47',
               'sele',
               'cdh5',
               'pecam1a',
               'mcamb',
               'podxl',
               'lef1',
               'fzd3b',
               'notum1a',
               'axin2',
               'LOC111580661',
               'tnfrsf19',
               'alpl',
               'LOC111576983',
               'cd36',
               'cd151',
               'LOC111582684',
               'esamb',
               'il13ra1',
               'itga4',
               'LOC111584896',
               'LOC111574700',
               'ralbp1'))+
  coord_flip()


DotPlot(obj, c(   'LOC111588076', #gad1,
                'gad2',
                'LOC111577366',
               'pvalb3',
               'pvalb5',
               'pvalb6',
               'pvalb7',#i also want to look at interneurons
               'pvalb8',
               'pvalb9',
               'LOC111575367',
               'LOC111574793', #C1ql1
               'dock10',
               'gfra1b',
               'plch2b',
               'ttll1',
               'tac1',
               'LOC111562384', #ccka
               'cckb',
               'sst1.1',
               'sst1.2'))+
  coord_flip()
