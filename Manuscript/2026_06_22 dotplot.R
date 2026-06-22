obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

genes = c(
  # cell types
  'elavl3',
  
  
  'slc17a6b',
  'slc17a7a', #glut
  'LOC111588076', #gad1b
  'gad2',#gaba 
  

  
  
  #0
  'soul5l', 
  'slc32a1',
  
  #1
  'LOC111577263',
  'gfap',
  
  #2
  'mbpa',
  
  #3
  'kcnt1a',
  'scn5lab',
  
  #4
  'tac1',
  'arpp21',
  
  #5
  'baiap3',
  'sema3h',
  
  #POA
  'prkd1',
  'hmx2',
  'hmx3a',
  'hmx3b',
  'avp',
  'oxt',
  'cckb', 
  'esr2b',
  'ar',
  
  #7
  'LOC111562384',
  'znf536',
  #8
  'LOC111568697',
  'LOC111575520',
  
  #9
  'mex3a',
  #10
  'LOC111582900',
  'prkcq',
  
  #11
  'ptprc',
  'p2ry12',
  
  #12
  'susd5',
  'ebf3b',
  
  #13
  'cspg4',
  
  # vsst
  'npy',
  'sst1.1',
  
  #15
  'crocc2',
  
  #16
  'rap1gapl',
  'nptx2b',
  
  #17
  'nr5a2',
  'calb2a',
  
  #18
  'LOC111583752',
  'ebf1a',
  
  #19
  'gata3',
  'gata2b',
  
  #21
  'lhx6a',
  'adarb2',
  
  #22
  'LOC111583752',
  'zic2a',
  
  #23
  'pde9ac',
  'chgb',
  
  #24
  'fat2',
  'gabra6b',
  
  #25
  'meis2a',
  'cart2',
  
  #26
  'cenpp',
  'top2a'
)%>%unique()

DotPlot(obj, genes, dot.min = 0.1)+
  coord_flip()
