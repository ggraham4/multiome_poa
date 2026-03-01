phase_degs = read.csv('~/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

phase_degs$d_m_sig = ifelse(phase_degs$d_m_p.value<0.05, '*','Not')
phase_degs$d_f_sig = ifelse(phase_degs$d_f_p.value<0.05, '*','Not')
phase_degs$f_m_sig = ifelse(phase_degs$f_m_p.value<0.05, '*','Not')

table(phase_degs$d_m_sig)
table(phase_degs$d_f_sig)
table(phase_degs$f_m_sig)

phase_degs_10=subset(phase_degs, cluster ==10)


time_dict = c(
  'Early' = 1,
  'Transiently' = 2,
  'Late' = 3,
  'Progressively' =2
)

phase_degs$time_number = time_dict[phase_degs$first_word]

phase_degs_timed =phase_degs%>%
  group_by(cluster)%>%
  summarize(time = mean(time_number),
            n_degs = n())
# cluster 10 is one of the first major ones, then 5, 2, 11

ggplot(phase_degs_timed, aes(x = time, y = n_degs))+
  geom_text(aes(label = cluster))

DotPlot(obj, c('LOC111571064',
               'gnrh2',
               'gnrh3',
               'gnrhr1',
               'gnrhr4'))+
  coord_flip()

DotPlot(obj, c('otpa',
               'fezf2',
               'LOC111564688', #otpb-like
               'slc6a3',
               'th',
               'th2',
               'sst1.1', #sst1.1+
               'avp',
               'oxt',
               'LOC111571064'))+
  coord_flip()
#dopamine 


DotPlot(obj, c('crhb',
               'isl1a',
               'dlx5a',
               'drd3',
               'LOC111566447',
               'LOC111575776', # penka
               'cckb',
               'LOC111562384',#ccka
               'vip',
               'trh',
               'nts',
               'oxt',
               'avp',
               'otpa',
               'fezf2',
               'arnt2','nkx2.1','shha',
               'tbx2b',
               'rx3',
               'gal',
               'kiss1',
               'hmx2',
               'hmx3a',
               'prok2',
               'npy',
               'pmch',
               'lhx2b',
               'six6a',
               'LOC111584656' # pomc

               ))+
  coord_flip()

