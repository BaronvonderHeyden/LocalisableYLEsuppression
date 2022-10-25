cd("C:\\Users\\rene4\\OneDrive\\Research Project London\\GitHub")

sf72 = readdlm("sf_scan72", ',')
sf73 = readdlm("sf_scan73", ',')
sf_lowRm78 =  readdlm("sf_scan78_veryLowRmOf2", ',')
sf_lowRm79 =  readdlm("sf_scan79_veryLowRmOf2", ',')
sf_highRm710 =  readdlm("sf_scan710_veryHighRmOf12", ',')
sf_highRm711  =  readdlm("sf_scan711_veryHighRmOf12", ',')

sf_scan = sensitivity_plot3(sf72, sf73, sf_lowRm78, sf_lowRm79, sf_highRm710, sf_highRm711)
savefig(sf_scan, "sf_scan_extreme.svg")


#
# sf_lowRm78 =  readdlm("sf_scan74", ',')
# sf_lowRm79 =  readdlm("sf_scan75", ',')
# sf_highRm710 =  readdlm("sf_scan76", ',')
# sf_highRm711  =  readdlm("sf_scan77", ',')
#
#
# sf_scan = sensitivity_plot3(sf72, sf73, sf_lowRm78, sf_lowRm79, sf_highRm710, sf_highRm711)
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(sf_scan, "sf_scan_moderate.svg")

#---------------------------

hf72 = readdlm("hf_scan72", ',')
hf73 = readdlm("hf_scan73", ',')
hf_lowRm78 =  readdlm("hf_scan78_veryLowRmOf2", ',')
hf_lowRm79 =  readdlm("hf_scan79_veryLowRmOf2", ',')
hf_highRm710 =  readdlm("hf_scan710_veryHighRmOf12", ',')
hf_highRm711  =  readdlm("hf_scan711_veryHighRmOf12", ',')

hf_scan = sensitivity_plot3(hf72, hf73, hf_lowRm78, hf_lowRm79, hf_highRm710, hf_highRm711)
savefig(hf_scan, "hf_scan_extreme.svg")


#fat grey line is below ymin
# plot(hf72[2:end,1], (hf_lowRm78[2:end,2]), colour=:grey, dpi = 300, legend = false, yaxis =:log, ylims = (1e-40,20))


# hf_lowRm78 =  readdlm("hf_scan74", ',')
# hf_lowRm79 =  readdlm("hf_scan75", ',')
# hf_highRm710 =  readdlm("hf_scan76", ',')
# hf_highRm711  =  readdlm("hf_scan77", ',')
#
#
# hf_scan = sensitivity_plot3(hf72, hf73, hf_lowRm78, hf_lowRm79, hf_highRm710, hf_highRm711)
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(hf_scan, "hf_scan_moderate.svg")


#-----------------------------------------


sa72 = readdlm("sa_scan72", ',')
sa73 = readdlm("sa_scan73", ',')

sa_lowRm78 =  readdlm("sa_scan78_veryLowRmOf2", ',')
sa_lowRm79 =  readdlm("sa_scan79_veryLowRmOf2", ',')
sa_highRm710 =  readdlm("sa_scan710_veryHighRmOf12", ',')
sa_highRm711  =  readdlm("sa_scan711_veryHighRmOf12", ',')

sa_scan = sensitivity_plot3(sa72, sa73, sa_lowRm78, sa_lowRm79, sa_highRm710, sa_highRm711)
plot!(xlims=(0,0.5))
savefig(sa_scan, "sa_scan_extreme.svg")


#
# sa_lowRm78 =  readdlm("sa_scan74", ',')
# sa_lowRm79 =  readdlm("sa_scan75", ',')
# sa_highRm710 =  readdlm("sa_scan76", ',')
# sa_highRm711  =  readdlm("sa_scan77", ',')
#
#
# sa_scan = sensitivity_plot3(sa72, sa73, sa_lowRm78, sa_lowRm79, sa_highRm710, sa_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(sa_scan, "sa_scan_moderate.svg")








#-------------------------------------------------------------

sb72 = readdlm("sb_scan72", ',')
sb73 = readdlm("sb_scan73", ',')


sb_lowRm78 =  readdlm("sb_scan78_veryLowRmOf2", ',')
sb_lowRm79 =  readdlm("sb_scan79_veryLowRmOf2", ',')
sb_highRm710 =  readdlm("sb_scan710_veryHighRmOf12", ',')
sb_highRm711  =  readdlm("sb_scan711_veryHighRmOf12", ',')

sb_scan = sensitivity_plot3(sb72, sb73, sb_lowRm78, sb_lowRm79, sb_highRm710, sb_highRm711)
plot!(xlims=(0,0.5))
savefig(sb_scan, "sb_scan_extreme.svg")

#
#
# sb_lowRm78 =  readdlm("sb_scan74", ',')
# sb_lowRm79 =  readdlm("sb_scan75", ',')
# sb_highRm710 =  readdlm("sb_scan76", ',')
# sb_highRm711  =  readdlm("sb_scan77", ',')
#
#
# sb_scan = sensitivity_plot3(sb72, sb73, sb_lowRm78, sb_lowRm79, sb_highRm710, sb_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(sb_scan, "sb_scan_moderate.svg")
#




#----------------------------------------------------


sc72 = readdlm("sc_scan72", ',')
sc73 = readdlm("sc_scan73", ',')


sc_lowRm78 =  readdlm("sc_scan78_veryLowRmOf2", ',')
sc_lowRm79 =  readdlm("sc_scan79_veryLowRmOf2", ',')
sc_highRm710 =  readdlm("sc_scan710_veryHighRmOf12", ',')
sc_highRm711  =  readdlm("sc_scan711_veryHighRmOf12", ',')

sc_scan = sensitivity_plot3(sc72, sc73, sc_lowRm78, sc_lowRm79, sc_highRm710, sc_highRm711)
plot!(xlims=(0,0.5))
savefig(sc_scan, "sc_scan_extreme.svg")

#
#
# sc_lowRm78 =  readdlm("sc_scan74", ',')
# sc_lowRm79 =  readdlm("sc_scan75", ',')
# sc_highRm710 =  readdlm("sc_scan76", ',')
# sc_highRm711  =  readdlm("sc_scan77", ',')
#
#
# sc_scan = sensitivity_plot3(sc72, sc73, sc_lowRm78, sc_lowRm79, sc_highRm710, sc_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(sc_scan, "sc_scan_moderate.svg")







#-------------------------------------------------------

sd72 = readdlm("sd_scan72", ',')
sd73 = readdlm("sd_scan73", ',')


sd_lowRm78 =  readdlm("sd_scan78_veryLowRmOf2", ',')
sd_lowRm79 =  readdlm("sd_scan79_veryLowRmOf2", ',')
sd_highRm710 =  readdlm("sd_scan710_veryHighRmOf12", ',')
sd_highRm711  =  readdlm("sd_scan711_veryHighRmOf12", ',')

sd_scan = sensitivity_plot3(sd72, sd73, sd_lowRm78, sd_lowRm79, sd_highRm710, sd_highRm711)
plot!(xlims=(0,0.5))
savefig(sd_scan, "sd_scan_extreme.svg")


#
# sd_lowRm78 =  readdlm("sd_scan74", ',')
# sd_lowRm79 =  readdlm("sd_scan75", ',')
# sd_highRm710 =  readdlm("sd_scan76", ',')
# sd_highRm711  =  readdlm("sd_scan77", ',')
#
#
# sd_scan = sensitivity_plot3(sd72, sd73, sd_lowRm78, sd_lowRm79, sd_highRm710, sd_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(sd_scan, "sd_scan_moderate.svg")
#







#-----------------------------------------------------

se72 = readdlm("se_scan72", ',')
se73 = readdlm("se_scan73", ',')


se_lowRm78 =  readdlm("se_scan78_veryLowRmOf2", ',')
se_lowRm79 =  readdlm("se_scan79_veryLowRmOf2", ',')
se_highRm710 =  readdlm("se_scan710_veryHighRmOf12", ',')
se_highRm711  =  readdlm("se_scan711_veryHighRmOf12", ',')

se_scan = sensitivity_plot3(se72, se73, se_lowRm78, se_lowRm79, se_highRm710, se_highRm711)
plot!(xlims=(0,0.5))
savefig(se_scan, "se_scan_extreme.svg")

#
#
# se_lowRm78 =  readdlm("se_scan74", ',')
# se_lowRm79 =  readdlm("se_scan75", ',')
# se_highRm710 =  readdlm("se_scan76", ',')
# se_highRm711  =  readdlm("se_scan77", ',')
#
#
# se_scan = sensitivity_plot3(se72, se73, se_lowRm78, se_lowRm79, se_highRm710, se_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(se_scan, "se_scan_moderate.svg")





#--------------------------------------------------------
#
# he72 = readdlm("he_scan72", ',')
# he73 = readdlm("he_scan73", ',')
#
# he_lowRm78 =  readdlm("he_scan78_veryLowRmOf2", ',')
# he_lowRm79 =  readdlm("he_scan79_veryLowRmOf2", ',')
# he_highRm710 =  readdlm("he_scan710_veryHighRmOf12", ',')
# he_highRm711  =  readdlm("he_scan711_veryHighRmOf12", ',')
#
# he_scan = sensitivity_plot3(he72, he73, he_lowRm78, he_lowRm79, he_highRm710, he_highRm711)
# plot!(xlims=(0,0.5))
# savefig(he_scan, "he_scan_extreme.svg")
#
#
#
# he_lowRm78 =  readdlm("he_scan74", ',')
# he_lowRm79 =  readdlm("he_scan75", ',')
# he_highRm710 =  readdlm("he_scan76", ',')
# he_highRm711  =  readdlm("he_scan77", ',')
#
#
# he_scan = sensitivity_plot3(he72, he73, he_lowRm78, he_lowRm79, he_highRm710, he_highRm711)
# plot!(xlims=(0,0.5))
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(he_scan, "he_scan_moderate.svg")
#






#------------------------------------------------------

#
# sm72 = readdlm("sm_scan72", ',')
# sm73 = readdlm("sm_scan73", ',')


ee72 = readdlm("ee_scan72", ',')
ee72[2,2] = NaN
ee73 = readdlm("ee_scan73", ',')


ee_lowRm78 =  readdlm("ee_scan78_veryLowRmOf2", ',')
ee_lowRm79 =  readdlm("ee_scan79_veryLowRmOf2", ',')
ee_highRm710 =  readdlm("ee_scan710_veryHighRmOf12", ',')
ee_highRm711  =  readdlm("ee_scan711_veryHighRmOf12", ',')

ee_lowRm78[2,2] = NaN
ee_lowRm79[2,2] = NaN
ee_highRm710[2,2] = NaN

ee_scan = sensitivity_plot3(ee72, ee73, ee_lowRm78, ee_lowRm79, ee_highRm710, ee_highRm711)

savefig(ee_scan, "ee_scan_extreme.svg")

#
#
# ee_lowRm78 =  readdlm("ee_scan74", ',')
# ee_lowRm79 =  readdlm("ee_scan75", ',')
# ee_highRm710 =  readdlm("ee_scan76", ',')
# ee_highRm711  =  readdlm("ee_scan77", ',')
#
# ee_lowRm78[2,2] = NaN
# ee_lowRm79[2,2] = NaN
# ee_highRm710[2,2] = NaN
#
#
# ee_scan = sensitivity_plot3(ee72, ee73, ee_lowRm78, ee_lowRm79, ee_highRm710, ee_highRm711)
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(ee_scan, "ee_scan_moderate.svg")


#-----------------------------------------------


es72 = readdlm("es_scan72", ',')
es73 = readdlm("es_scan73", ',')

es_lowRm78 =  readdlm("es_scan78_veryLowRmOf2", ',')
es_lowRm79 =  readdlm("es_scan79_veryLowRmOf2", ',')
es_highRm710 =  readdlm("es_scan710_veryHighRmOf12", ',')
es_highRm711  =  readdlm("es_scan711_veryHighRmOf12", ',')

es_scan = sensitivity_plot3(es72, es73, es_lowRm78, es_lowRm79, es_highRm710, es_highRm711)

savefig(es_scan, "es_scan_extreme.svg")

#
#
# es_lowRm78 =  readdlm("es_scan74", ',')
# es_lowRm79 =  readdlm("es_scan75", ',')
# es_highRm710 =  readdlm("es_scan76", ',')
# es_highRm711  =  readdlm("es_scan77", ',')
#
#
# es_scan = sensitivity_plot3(es72, es73, es_lowRm78, es_lowRm79, es_highRm710, es_highRm711)
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(es_scan, "es_scan_moderate.svg")
#
#



#------------------------------------------------

eh72 = readdlm("eh_scan72", ',')
eh73 = readdlm("eh_scan73", ',')


eh_lowRm78 =  readdlm("eh_scan78_veryLowRmOf2", ',')
eh_lowRm79 =  readdlm("eh_scan79_veryLowRmOf2", ',')
eh_highRm710 =  readdlm("eh_scan710_veryHighRmOf12", ',')
eh_highRm711  =  readdlm("eh_scan711_veryHighRmOf12", ',')

eh_scan = sensitivity_plot3(eh72, eh73, eh_lowRm78, eh_lowRm79, eh_highRm710, eh_highRm711)
savefig(eh_scan, "eh_scan_extreme.svg")


#
# eh_lowRm78 =  readdlm("eh_scan74", ',')
# eh_lowRm79 =  readdlm("eh_scan75", ',')
# eh_highRm710 =  readdlm("eh_scan76", ',')
# eh_highRm711  =  readdlm("eh_scan77", ',')
#
#
# eh_scan = sensitivity_plot3(eh72, eh73, eh_lowRm78, eh_lowRm79, eh_highRm710, eh_highRm711)
# plot!(ytickfontsize = 14, xtickfontsize = 14)
# savefig(eh_scan, "eh_scan_moderate.svg")
#




#------------------------------------------------


# he272 = readdlm("he2_scan72", ',')
# he273 = readdlm("he2_scan73", ',')


m172 = readdlm("m1_scan72", ',')
m173 = readdlm("m1_scan73", ',')

m1_lowRm78 =  readdlm("m1_scan78_veryLowRmOf2", ',')
m1_lowRm79 =  readdlm("m1_scan79_veryLowRmOf2", ',')
m1_highRm710 =  readdlm("m1_scan710_veryHighRmOf12", ',')
m1_highRm711  =  readdlm("m1_scan711_veryHighRmOf12", ',')

# m1 values are stored as exponent to base 10 values, so we first transform values

for i in 2:length(m172[:,1])
    m172[i,1] = 10^(m172[i,1])
    m173[i,1] = 10^(m173[i,1])
    m1_lowRm78[i,1] = 10^(m1_lowRm78[i,1])
    m1_lowRm79[i,1] = 10^(m1_lowRm79[i,1])
    m1_highRm710[i,1] = 10^(m1_highRm710[i,1])
    m1_highRm711[i,1]  = 10^(m1_highRm711[i,1])

end

m1_scan = sensitivity_plot3(m172, m173, m1_lowRm78, m1_lowRm79, m1_highRm710, m1_highRm711)
plot!(xaxis=:log, xlims = (1e-3,1), ytickfontsize = 14, xtickfontsize = 14)
savefig(m1_scan, "m1_scan_extreme.svg")

#
#
# m1_lowRm78 =  readdlm("m1_scan74", ',')
# m1_lowRm79 =  readdlm("m1_scan75", ',')
# m1_highRm710 =  readdlm("m1_scan76", ',')
# m1_highRm711  =  readdlm("m1_scan77", ',')
#
# # m1 values are stored as exponent to base 10 values, so we first transform values
#
# for i in 2:length(m172[:,1])
#     m1_lowRm78[i,1] = 10^(m1_lowRm78[i,1])
#     m1_lowRm79[i,1] = 10^(m1_lowRm79[i,1])
#     m1_highRm710[i,1] = 10^(m1_highRm710[i,1])
#     m1_highRm711[i,1]  = 10^(m1_highRm711[i,1])
#
# end
#
# m1_scan = sensitivity_plot3(m172, m173, m1_lowRm78, m1_lowRm79, m1_highRm710, m1_highRm711)
# plot!(xaxis=:log, xlims = (1e-3,1), ytickfontsize = 14, xtickfontsize = 14)
# savefig(m1_scan, "m1_scan_moderate.svg")
#
#
#





#-------------------------------------------------


m272 = readdlm("m2_scan72", ',')
m273 = readdlm("m2_scan73", ',')


m2_lowRm78 =  readdlm("m2_scan78_veryLowRmOf2", ',')
m2_lowRm79 =  readdlm("m2_scan79_veryLowRmOf2", ',')
m2_highRm710 =  readdlm("m2_scan710_veryHighRmOf12", ',')
m2_highRm711  =  readdlm("m2_scan711_veryHighRmOf12", ',')

# m1 values are stored as exponent to base 10 values, so we first transform values

for i in 2:length(m272[:,1])
    m272[i,1] = 10^(m272[i,1])
    m273[i,1] = 10^(m273[i,1])
    m2_lowRm78[i,1] = 10^(m2_lowRm78[i,1])
    m2_lowRm79[i,1] = 10^(m2_lowRm79[i,1])
    m2_highRm710[i,1] = 10^(m2_highRm710[i,1])
    m2_highRm711[i,1]  = 10^(m2_highRm711[i,1])

end

m2_scan = sensitivity_plot3(m272, m273, m2_lowRm78, m2_lowRm79, m2_highRm710, m2_highRm711)
plot!(xaxis=:log, xlims = (1e-3,1), ytickfontsize = 14, xtickfontsize = 14)
savefig(m2_scan, "m2_scan_extreme.svg")

#
#
# m2_lowRm78 =  readdlm("m2_scan74", ',')
# m2_lowRm79 =  readdlm("m2_scan75", ',')
# m2_highRm710 =  readdlm("m2_scan76", ',')
# m2_highRm711  =  readdlm("m2_scan77", ',')
#
# # m1 values are stored as exponent to base 10 values, so we first transform values
#
# for i in 2:length(m272[:,1])
#     m2_lowRm78[i,1] = 10^(m2_lowRm78[i,1])
#     m2_lowRm79[i,1] = 10^(m2_lowRm79[i,1])
#     m2_highRm710[i,1] = 10^(m2_highRm710[i,1])
#     m2_highRm711[i,1]  = 10^(m2_highRm711[i,1])
#
# end
#
# m2_scan = sensitivity_plot3(m272, m273, m2_lowRm78, m2_lowRm79, m2_highRm710, m2_highRm711)
# plot!(xaxis=:log, xlims = (1e-3,1), ytickfontsize = 14, xtickfontsize = 14)
# savefig(m2_scan, "m2_scan_moderate.svg")
#
#




#-------------------------------------------


er172 = readdlm("er1_scan72", ',')
er173 = readdlm("er1_scan73", ',')


er272 = readdlm("er2_scan72", ',')
er273 = readdlm("er2_scan73", ',')


er372 = readdlm("er3_scan72", ',')
er373 = readdlm("er3_scan73", ',')

Rm72 = readdlm("Rm_scan72", ',')
Rm73 = readdlm("Rm_scan73", ',')
