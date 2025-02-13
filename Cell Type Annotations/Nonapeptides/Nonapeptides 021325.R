DotPlot(object = obj,
        features= c('elavl3','oxt','avp'),
        dot.min =0.24,
        dot.scale=T)+
  coord_flip()

library(scCustomize)

avp_per = t(Percent_Expressing(obj, 'avp'))[,1]
names(avp_per) = 0:31
sort(avp_per)

avp_expression = AverageExpression(obj, features = 'avp', slot = 'data')[[1]]@x
names(avp_expression) = 0:31
sort(avp_expression)


oxt_per = t(Percent_Expressing(obj, 'oxt'))[,1]
names(oxt_per) = 0:31
sort(oxt_per)

oxt_expression = AverageExpression(obj, features = 'oxt', slot = 'data')[[1]]@x
names(oxt_expression) = 0:31
sort(oxt_expression)

