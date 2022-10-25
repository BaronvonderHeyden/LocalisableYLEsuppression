cd("C:/Users/rene4/OneDrive/Research Project London/GitHub")


#another parameter set
Parameters_set_sensitivity8_verylowRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("2"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.05"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity9_verylowRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("2"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.6"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#another parameter set
Parameters_set_sensitivity10_veryhighRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("12"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.05"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity11_veryhighRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("12"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.6"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)




sf = Any["value" "sup" "dur" "Y" "S"]
hf = Any["value" "sup" "dur" "Y" "S"]

sa = Any["value" "sup" "dur" "Y" "S"]
sb = Any["value" "sup" "dur" "Y" "S"]
sc = Any["value" "sup" "dur" "Y" "S"]
sd = Any["value" "sup" "dur" "Y" "S"]
se = Any["value" "sup" "dur" "Y" "S"]
he = Any["value" "sup" "dur" "Y" "S"]

prot_expression_cost = Any["value" "sup" "dur" "Y" "S"]

sm = Any["value" "sup" "dur" "Y" "S"]

ee = Any["value" "sup" "dur" "Y" "S"]
eh = Any["value" "sup" "dur" "Y" "S"]
es = Any["value" "sup" "dur" "Y" "S"]
he2 = Any["value" "sup" "dur" "Y" "S"]

shredding_activity_cost = Any["value" "sup" "dur" "Y" "S"]

m1 = Any["value" "sup" "dur" "Y" "S"]
m2 = Any["value" "sup" "dur" "Y" "S"]

er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]




apply_parameters_set(Parameters_set_sensitivity8_verylowRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity8_verylowRM)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global sf
    current_Parameters["s_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sf = vcat(sf, add)

end

writedlm("sf_scan78_veryLowRmOf2",  sf, ',')

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global hf
    current_Parameters["h_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    hf = vcat(hf, add)

end
writedlm("hf_scan78_veryLowRmOf2",  hf, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sa
    current_Parameters["s_a"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sa = vcat(sa, add)

end
writedlm("sa_scan78_veryLowRmOf2",  sa, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sb
    current_Parameters["s_b"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sb = vcat(sb, add)

end
writedlm("sb_scan78_veryLowRmOf2",  sb, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global sc
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sc = vcat(sc, add)

end
writedlm("sc_scan78_veryLowRmOf2",  sc, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sd
    current_Parameters["s_d"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sd = vcat(sd, add)

end
writedlm("sd_scan78_veryLowRmOf2",  sd, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global se
    current_Parameters["s_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    se = vcat(se, add)


end
writedlm("se_scan78_veryLowRmOf2",  se, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global he
    current_Parameters["h_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he = vcat(he, add)

end
writedlm("he_scan78_veryLowRmOf2",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


#
# for i in 0:0.005:1
#     global prot_expression_cost
#     current_Parameters["s_a"] = i
#     current_Parameters["s_c"] = i
#     current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     prot_expression_cost = vcat(prot_expression_cost, add)
#
# end
# writedlm("prot_expression_cost_scan78_veryLowRmOf2",  prot_expression_cost, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.02:1
    global sm
    current_Parameters["s_m"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sm = vcat(sm, add)


end
writedlm("sm_scan78_veryLowRmOf2",  sm, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global ee
    current_Parameters["e_e"] = i

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    ee = vcat(ee, add)

end
writedlm("ee_scan78_veryLowRmOf2",  ee, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global eh
    current_Parameters["e_h"] = i

    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    eh = vcat(eh, add)
end
writedlm("eh_scan78_veryLowRmOf2",  eh, ',')




current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in 1:-0.01:0
    global es
    current_Parameters["e_s"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    es = vcat(es, add)

end
writedlm("es_scan78_veryLowRmOf2",  es, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global he2
    current_Parameters["h_e2"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he2 = vcat(he2, add)
end

writedlm("he2_scan78_veryLowRmOf2",  he2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global shredding_activity_cost
    current_Parameters["h_e"] = i
    current_Parameters["h_e2"] = i

    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    shredding_activity_cost = vcat(shredding_activity_cost, add)
end

writedlm("shredding_activity_cost_scan78_veryLowRmOf2",  shredding_activity_cost, ',')


current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -3:0.05:0
    global m1
    current_Parameters["m_1"] = 10^(i)
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m1 = vcat(m1, add)

end

writedlm("m1_scan78_veryLowRmOf2",  m1, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in -3:0.05:0
    global m2
    current_Parameters["m_2"] = 10^(i)
    mutation_matrix0 = deepcopy(mutation_matrix_num)

    for i in 1:length(mutation_matrix[:,1])
        for j in mutation_symbols[i]
            mutation_matrix0[i,j] = BigFloat(float.(mutation_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["mutation_matrix"] = mutation_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m2 = vcat(m2, add)

end

writedlm("m2_scan78_veryLowRmOf2",  m2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -6:0.05:0
    global er1
    current_Parameters["er_1"] = 10^(i)

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er1 = vcat(er1, add)

end

writedlm("er1_scan78_veryLowRmOf2",  er1, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er2
    current_Parameters["er_2"] = i
    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end
    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er2 = vcat(er2, add)
end
writedlm("er2_scan78_veryLowRmOf2",  er2, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er3
    current_Parameters["er_3"] = i
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er3 = vcat(er3, add)
end
writedlm("er3_scan78_veryLowRmOf2",  er3, ',')



#- -------


sf = Any["value" "sup" "dur" "Y" "S"]
hf = Any["value" "sup" "dur" "Y" "S"]

sa = Any["value" "sup" "dur" "Y" "S"]
sb = Any["value" "sup" "dur" "Y" "S"]
sc = Any["value" "sup" "dur" "Y" "S"]
sd = Any["value" "sup" "dur" "Y" "S"]
se = Any["value" "sup" "dur" "Y" "S"]
he = Any["value" "sup" "dur" "Y" "S"]

prot_expression_cost = Any["value" "sup" "dur" "Y" "S"]

sm = Any["value" "sup" "dur" "Y" "S"]

ee = Any["value" "sup" "dur" "Y" "S"]
eh = Any["value" "sup" "dur" "Y" "S"]
es = Any["value" "sup" "dur" "Y" "S"]
he2 = Any["value" "sup" "dur" "Y" "S"]

shredding_activity_cost = Any["value" "sup" "dur" "Y" "S"]

m1 = Any["value" "sup" "dur" "Y" "S"]
m2 = Any["value" "sup" "dur" "Y" "S"]

er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]




apply_parameters_set(Parameters_set_sensitivity9_verylowRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity9_verylowRM)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global sf
    current_Parameters["s_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sf = vcat(sf, add)

end

writedlm("sf_scan79_veryLowRmOf2",  sf, ',')

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global hf
    current_Parameters["h_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    hf = vcat(hf, add)

end
writedlm("hf_scan79_veryLowRmOf2",  hf, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sa
    current_Parameters["s_a"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sa = vcat(sa, add)

end
writedlm("sa_scan79_veryLowRmOf2",  sa, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sb
    current_Parameters["s_b"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sb = vcat(sb, add)

end
writedlm("sb_scan79_veryLowRmOf2",  sb, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global sc
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sc = vcat(sc, add)

end
writedlm("sc_scan79_veryLowRmOf2",  sc, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sd
    current_Parameters["s_d"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sd = vcat(sd, add)

end
writedlm("sd_scan79_veryLowRmOf2",  sd, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global se
    current_Parameters["s_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    se = vcat(se, add)


end
writedlm("se_scan79_veryLowRmOf2",  se, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global he
    current_Parameters["h_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he = vcat(he, add)

end
writedlm("he_scan79_veryLowRmOf2",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


#
# for i in 0:0.005:1
#     global prot_expression_cost
#     current_Parameters["s_a"] = i
#     current_Parameters["s_c"] = i
#     current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     prot_expression_cost = vcat(prot_expression_cost, add)
#
# end
# writedlm("prot_expression_cost_scan79_veryLowRmOf2",  prot_expression_cost, ',')
#
#

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.02:1
    global sm
    current_Parameters["s_m"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sm = vcat(sm, add)


end
writedlm("sm_scan79_veryLowRmOf2",  sm, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global ee
    current_Parameters["e_e"] = i

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    ee = vcat(ee, add)

end
writedlm("ee_scan79_veryLowRmOf2",  ee, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global eh
    current_Parameters["e_h"] = i

    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    eh = vcat(eh, add)
end
writedlm("eh_scan79_veryLowRmOf2",  eh, ',')




current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in 1:-0.01:0
    global es
    current_Parameters["e_s"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    es = vcat(es, add)

end
writedlm("es_scan79_veryLowRmOf2",  es, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global he2
    current_Parameters["h_e2"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he2 = vcat(he2, add)
end

writedlm("he2_scan79_veryLowRmOf2",  he2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global shredding_activity_cost
    current_Parameters["h_e"] = i
    current_Parameters["h_e2"] = i

    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    shredding_activity_cost = vcat(shredding_activity_cost, add)
end

writedlm("shredding_activity_cost_scan79_veryLowRmOf2",  shredding_activity_cost, ',')


current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -3:0.05:0
    global m1
    current_Parameters["m_1"] = 10^(i)
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m1 = vcat(m1, add)

end

writedlm("m1_scan79_veryLowRmOf2",  m1, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in -3:0.05:0
    global m2
    current_Parameters["m_2"] = 10^(i)
    mutation_matrix0 = deepcopy(mutation_matrix_num)

    for i in 1:length(mutation_matrix[:,1])
        for j in mutation_symbols[i]
            mutation_matrix0[i,j] = BigFloat(float.(mutation_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["mutation_matrix"] = mutation_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m2 = vcat(m2, add)

end

writedlm("m2_scan79_veryLowRmOf2",  m2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -6:0.05:0
    global er1
    current_Parameters["er_1"] = 10^(i)

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er1 = vcat(er1, add)

end

writedlm("er1_scan79_veryLowRmOf2",  er1, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er2
    current_Parameters["er_2"] = i
    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end
    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er2 = vcat(er2, add)
end
writedlm("er2_scan79_veryLowRmOf2",  er2, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0.5:0.01:1
    global er3
    current_Parameters["er_3"] = i
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er3 = vcat(er3, add)
end
writedlm("er3_scan79_veryLowRmOf2",  er3, ',')



# ------------------
# also Rm

# apply_parameters_set(Parameters_set_sensitivity4_lowRM)
#
# Rm = Any["value" "sup" "dur" "Y" "S"]
#
# for i in 3:0.5:20
#     global Rm
#     current_Parameters["Rm"] = i
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     Rm = vcat(Rm, add)
# end
#
# writedlm("Rm_scan78_veryLowRmOf2",  Rm, ',')
#


# -------------------
# apply_parameters_set(Parameters_set_sensitivity9_verylowRM)
#
# Rm = Any["value" "sup" "dur" "Y" "S"]
#
# for i in 3:0.5:20
#     global Rm
#     current_Parameters["Rm"] = i
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     Rm = vcat(Rm, add)
# end
#
# writedlm("Rm_scan79_veryLowRmOf2",  Rm, ',')




sf = Any["value" "sup" "dur" "Y" "S"]
hf = Any["value" "sup" "dur" "Y" "S"]

sa = Any["value" "sup" "dur" "Y" "S"]
sb = Any["value" "sup" "dur" "Y" "S"]
sc = Any["value" "sup" "dur" "Y" "S"]
sd = Any["value" "sup" "dur" "Y" "S"]
se = Any["value" "sup" "dur" "Y" "S"]
he = Any["value" "sup" "dur" "Y" "S"]

prot_expression_cost = Any["value" "sup" "dur" "Y" "S"]

sm = Any["value" "sup" "dur" "Y" "S"]

ee = Any["value" "sup" "dur" "Y" "S"]
eh = Any["value" "sup" "dur" "Y" "S"]
es = Any["value" "sup" "dur" "Y" "S"]
he2 = Any["value" "sup" "dur" "Y" "S"]

shredding_activity_cost = Any["value" "sup" "dur" "Y" "S"]

m1 = Any["value" "sup" "dur" "Y" "S"]
m2 = Any["value" "sup" "dur" "Y" "S"]

er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]




apply_parameters_set(Parameters_set_sensitivity10_veryhighRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity10_veryhighRM)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global sf
    current_Parameters["s_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sf = vcat(sf, add)

end

writedlm("sf_scan710_veryHighRmOf12",  sf, ',')

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global hf
    current_Parameters["h_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    hf = vcat(hf, add)

end
writedlm("hf_scan710_veryHighRmOf12",  hf, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sa
    current_Parameters["s_a"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sa = vcat(sa, add)

end
writedlm("sa_scan710_veryHighRmOf12",  sa, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sb
    current_Parameters["s_b"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sb = vcat(sb, add)

end
writedlm("sb_scan710_veryHighRmOf12",  sb, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global sc
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sc = vcat(sc, add)

end
writedlm("sc_scan710_veryHighRmOf12",  sc, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sd
    current_Parameters["s_d"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sd = vcat(sd, add)

end
writedlm("sd_scan710_veryHighRmOf12",  sd, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global se
    current_Parameters["s_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    se = vcat(se, add)


end
writedlm("se_scan710_veryHighRmOf12",  se, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global he
    current_Parameters["h_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he = vcat(he, add)

end
writedlm("he_scan710_veryHighRmOf12",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


#
# for i in 0:0.005:1
#     global prot_expression_cost
#     current_Parameters["s_a"] = i
#     current_Parameters["s_c"] = i
#     current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     prot_expression_cost = vcat(prot_expression_cost, add)
#
# end
# writedlm("prot_expression_cost_scan710_veryHighRmOf12",  prot_expression_cost, ',')
#


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.02:1
    global sm
    current_Parameters["s_m"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sm = vcat(sm, add)


end
writedlm("sm_scan710_veryHighRmOf12",  sm, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global ee
    current_Parameters["e_e"] = i

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    ee = vcat(ee, add)

end
writedlm("ee_scan710_veryHighRmOf12",  ee, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global eh
    current_Parameters["e_h"] = i

    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    eh = vcat(eh, add)
end
writedlm("eh_scan710_veryHighRmOf12",  eh, ',')




current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in 1:-0.01:0
    global es
    current_Parameters["e_s"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    es = vcat(es, add)

end
writedlm("es_scan710_veryHighRmOf12",  es, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global he2
    current_Parameters["h_e2"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he2 = vcat(he2, add)
end

writedlm("he2_scan710_veryHighRmOf12",  he2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global shredding_activity_cost
    current_Parameters["h_e"] = i
    current_Parameters["h_e2"] = i

    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    shredding_activity_cost = vcat(shredding_activity_cost, add)
end

writedlm("shredding_activity_cost_scan710_veryHighRmOf12",  shredding_activity_cost, ',')


current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -3:0.05:0
    global m1
    current_Parameters["m_1"] = 10^(i)
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m1 = vcat(m1, add)

end

writedlm("m1_scan710_veryHighRmOf12",  m1, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in -3:0.05:0
    global m2
    current_Parameters["m_2"] = 10^(i)
    mutation_matrix0 = deepcopy(mutation_matrix_num)

    for i in 1:length(mutation_matrix[:,1])
        for j in mutation_symbols[i]
            mutation_matrix0[i,j] = BigFloat(float.(mutation_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["mutation_matrix"] = mutation_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m2 = vcat(m2, add)

end

writedlm("m2_scan710_veryHighRmOf12",  m2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -6:0.05:0
    global er1
    current_Parameters["er_1"] = 10^(i)

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er1 = vcat(er1, add)

end

writedlm("er1_scan710_veryHighRmOf12",  er1, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er2
    current_Parameters["er_2"] = i
    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end
    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er2 = vcat(er2, add)
end
writedlm("er2_scan710_veryHighRmOf12",  er2, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er3
    current_Parameters["er_3"] = i
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er3 = vcat(er3, add)
end
writedlm("er3_scan710_veryHighRmOf12",  er3, ',')



#- -------


sf = Any["value" "sup" "dur" "Y" "S"]
hf = Any["value" "sup" "dur" "Y" "S"]

sa = Any["value" "sup" "dur" "Y" "S"]
sb = Any["value" "sup" "dur" "Y" "S"]
sc = Any["value" "sup" "dur" "Y" "S"]
sd = Any["value" "sup" "dur" "Y" "S"]
se = Any["value" "sup" "dur" "Y" "S"]
he = Any["value" "sup" "dur" "Y" "S"]

prot_expression_cost = Any["value" "sup" "dur" "Y" "S"]

sm = Any["value" "sup" "dur" "Y" "S"]

ee = Any["value" "sup" "dur" "Y" "S"]
eh = Any["value" "sup" "dur" "Y" "S"]
es = Any["value" "sup" "dur" "Y" "S"]
he2 = Any["value" "sup" "dur" "Y" "S"]

shredding_activity_cost = Any["value" "sup" "dur" "Y" "S"]

m1 = Any["value" "sup" "dur" "Y" "S"]
m2 = Any["value" "sup" "dur" "Y" "S"]

er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]




apply_parameters_set(Parameters_set_sensitivity11_veryhighRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity11_veryhighRM)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global sf
    current_Parameters["s_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sf = vcat(sf, add)

end

writedlm("sf_scan711_veryHighRmOf12",  sf, ',')

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global hf
    current_Parameters["h_f"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    hf = vcat(hf, add)

end
writedlm("hf_scan711_veryHighRmOf12",  hf, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sa
    current_Parameters["s_a"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sa = vcat(sa, add)

end
writedlm("sa_scan711_veryHighRmOf12",  sa, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sb
    current_Parameters["s_b"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sb = vcat(sb, add)

end
writedlm("sb_scan711_veryHighRmOf12",  sb, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global sc
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sc = vcat(sc, add)

end
writedlm("sc_scan711_veryHighRmOf12",  sc, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global sd
    current_Parameters["s_d"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sd = vcat(sd, add)

end
writedlm("sd_scan711_veryHighRmOf12",  sd, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global se
    current_Parameters["s_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    se = vcat(se, add)


end
writedlm("se_scan711_veryHighRmOf12",  se, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.005:1
    global he
    current_Parameters["h_e"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he = vcat(he, add)

end
writedlm("he_scan711_veryHighRmOf12",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


#
# for i in 0:0.005:1
#     global prot_expression_cost
#     current_Parameters["s_a"] = i
#     current_Parameters["s_c"] = i
#     current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     prot_expression_cost = vcat(prot_expression_cost, add)
#
# end
# writedlm("prot_expression_cost_scan711_veryHighRmOf12",  prot_expression_cost, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.02:1
    global sm
    current_Parameters["s_m"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    sm = vcat(sm, add)


end
writedlm("sm_scan711_veryHighRmOf12",  sm, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global ee
    current_Parameters["e_e"] = i

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    ee = vcat(ee, add)

end
writedlm("ee_scan711_veryHighRmOf12",  ee, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global eh
    current_Parameters["e_h"] = i

    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    eh = vcat(eh, add)
end
writedlm("eh_scan711_veryHighRmOf12",  eh, ',')




current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in 1:-0.01:0
    global es
    current_Parameters["e_s"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    es = vcat(es, add)

end
writedlm("es_scan711_veryHighRmOf12",  es, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global he2
    current_Parameters["h_e2"] = i

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    he2 = vcat(he2, add)
end

writedlm("he2_scan711_veryHighRmOf12",  he2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 1:-0.01:0
    global shredding_activity_cost
    current_Parameters["h_e"] = i
    current_Parameters["h_e2"] = i

    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    shredding_activity_cost = vcat(shredding_activity_cost, add)
end

writedlm("shredding_activity_cost_scan711_veryHighRmOf12",  shredding_activity_cost, ',')


current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -3:0.05:0
    global m1
    current_Parameters["m_1"] = 10^(i)
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m1 = vcat(m1, add)

end

writedlm("m1_scan711_veryHighRmOf12",  m1, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)

for i in -3:0.05:0
    global m2
    current_Parameters["m_2"] = 10^(i)
    mutation_matrix0 = deepcopy(mutation_matrix_num)

    for i in 1:length(mutation_matrix[:,1])
        for j in mutation_symbols[i]
            mutation_matrix0[i,j] = BigFloat(float.(mutation_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["mutation_matrix"] = mutation_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    m2 = vcat(m2, add)

end

writedlm("m2_scan711_veryHighRmOf12",  m2, ',')

current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in -6:0.05:0
    global er1
    current_Parameters["er_1"] = 10^(i)

    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["editing_matrix"] = editing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er1 = vcat(er1, add)

end

writedlm("er1_scan711_veryHighRmOf12",  er1, ',')



current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global er2
    current_Parameters["er_2"] = i
    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end
    current_matrices["sperm_matrix"] = sperm_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er2 = vcat(er2, add)
end
writedlm("er2_scan711_veryHighRmOf12",  er2, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0.5:0.01:1
    global er3
    current_Parameters["er_3"] = i
    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end

    current_matrices["homing_matrix"] = homing_matrix0
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    er3 = vcat(er3, add)
end
writedlm("er3_scan711_veryHighRmOf12",  er3, ',')



# ------------------
# also Rm
#
# apply_parameters_set(Parameters_set_sensitivity10_veryhighRM)
#
# Rm = Any["value" "sup" "dur" "Y" "S"]
#
# for i in 3:0.5:20
#     global Rm
#     current_Parameters["Rm"] = i
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     Rm = vcat(Rm, add)
# end
#
# writedlm("Rm_scan710_veryHighRmOf12",  Rm, ',')
#
#
#
# # -------------------
# apply_parameters_set(Parameters_set_sensitivity11_veryhighRM)
#
# Rm = Any["value" "sup" "dur" "Y" "S"]
#
# for i in 3:0.5:20
#     global Rm
#     current_Parameters["Rm"] = i
#     sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
#     add = [i sup dur Y S]
#     Rm = vcat(Rm, add)
# end
#
# writedlm("Rm_scan711_veryHighRmOf12",  Rm, ',')
