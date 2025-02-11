genes <- c('elavl3',
           'gad2', #GABA
           'LOC111588076', #GABA
           'LOC111584103', #GLU
           'slc17a6b', #glu
           'tac1',
           'tacr1a',
          'tacr1b',
          'tacr2',
           'tac3a',
          'tacr3a',
          'tacr3l',
          'npy',
          'kiss1',
          'kiss2',
          'kiss1ra',
          'kiss1rb',
          'esr1',
          'esr2a',
          'esr2b',
          'ar',
          'LOC111562384', #ccka
          'cckar',
          'cckb',
          'cckbra',
          'cckbrb',
          'gal',
          'galr1a',
          'galr1b',
          'LOC111566801',# galr2b
          'oxt',
          'avp'
          )

DotPlot(object = cluster_19, 
                 group.by = "sub", 
                 features = genes
        ) + 
  coord_flip()+
  geom_hline(yintercept = 20, linetype = 2, alpha = 0.5)
