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




apply_parameters_set(Parameters_set_sensitivity4_lowRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity4_lowRM)

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

writedlm("sf_scan74",  sf, ',')

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
writedlm("hf_scan74",  hf, ',')


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
writedlm("sa_scan74",  sa, ',')


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
writedlm("sb_scan74",  sb, ',')


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
writedlm("sc_scan74",  sc, ',')



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
writedlm("sd_scan74",  sd, ',')



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
writedlm("se_scan74",  se, ',')



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
writedlm("he_scan74",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global prot_expression_cost
    current_Parameters["s_a"] = i
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    prot_expression_cost = vcat(prot_expression_cost, add)

end
writedlm("prot_expression_cost_scan74",  prot_expression_cost, ',')



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
writedlm("sm_scan74",  sm, ',')


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
writedlm("ee_scan74",  ee, ',')



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
writedlm("eh_scan74",  eh, ',')




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
writedlm("es_scan74",  es, ',')


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

writedlm("he2_scan74",  he2, ',')

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

writedlm("shredding_activity_cost_scan74",  shredding_activity_cost, ',')


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

writedlm("m1_scan74",  m1, ',')

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

writedlm("m2_scan74",  m2, ',')

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

writedlm("er1_scan74",  er1, ',')



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
writedlm("er2_scan74",  er2, ',')


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
writedlm("er3_scan74",  er3, ',')



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




apply_parameters_set(Parameters_set_sensitivity5_lowRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity5_lowRM)

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

writedlm("sf_scan75",  sf, ',')

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
writedlm("hf_scan75",  hf, ',')


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
writedlm("sa_scan75",  sa, ',')


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
writedlm("sb_scan75",  sb, ',')


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
writedlm("sc_scan75",  sc, ',')



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
writedlm("sd_scan75",  sd, ',')



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
writedlm("se_scan75",  se, ',')



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
writedlm("he_scan75",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global prot_expression_cost
    current_Parameters["s_a"] = i
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    prot_expression_cost = vcat(prot_expression_cost, add)

end
writedlm("prot_expression_cost_scan75",  prot_expression_cost, ',')



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
writedlm("sm_scan75",  sm, ',')


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
writedlm("ee_scan75",  ee, ',')



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
writedlm("eh_scan75",  eh, ',')




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
writedlm("es_scan75",  es, ',')


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

writedlm("he2_scan75",  he2, ',')

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

writedlm("shredding_activity_cost_scan75",  shredding_activity_cost, ',')


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

writedlm("m1_scan75",  m1, ',')

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

writedlm("m2_scan75",  m2, ',')

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

writedlm("er1_scan75",  er1, ',')



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
writedlm("er2_scan75",  er2, ',')


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
writedlm("er3_scan75",  er3, ',')



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
# writedlm("Rm_scan74",  Rm, ',')
#


# -------------------
# apply_parameters_set(Parameters_set_sensitivity5_lowRM)
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
# writedlm("Rm_scan75",  Rm, ',')




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




apply_parameters_set(Parameters_set_sensitivity6_highRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity6_highRM)

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

writedlm("sf_scan76",  sf, ',')

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
writedlm("hf_scan76",  hf, ',')


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
writedlm("sa_scan76",  sa, ',')


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
writedlm("sb_scan76",  sb, ',')


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
writedlm("sc_scan76",  sc, ',')



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
writedlm("sd_scan76",  sd, ',')



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
writedlm("se_scan76",  se, ',')



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
writedlm("he_scan76",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global prot_expression_cost
    current_Parameters["s_a"] = i
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    prot_expression_cost = vcat(prot_expression_cost, add)

end
writedlm("prot_expression_cost_scan76",  prot_expression_cost, ',')



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
writedlm("sm_scan76",  sm, ',')


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
writedlm("ee_scan76",  ee, ',')



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
writedlm("eh_scan76",  eh, ',')




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
writedlm("es_scan76",  es, ',')


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

writedlm("he2_scan76",  he2, ',')

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

writedlm("shredding_activity_cost_scan76",  shredding_activity_cost, ',')


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

writedlm("m1_scan76",  m1, ',')

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

writedlm("m2_scan76",  m2, ',')

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

writedlm("er1_scan76",  er1, ',')



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
writedlm("er2_scan76",  er2, ',')


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
writedlm("er3_scan76",  er3, ',')



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




apply_parameters_set(Parameters_set_sensitivity7_highRM)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity7_highRM)

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

writedlm("sf_scan77",  sf, ',')

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
writedlm("hf_scan77",  hf, ',')


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
writedlm("sa_scan77",  sa, ',')


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
writedlm("sb_scan77",  sb, ',')


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
writedlm("sc_scan77",  sc, ',')



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
writedlm("sd_scan77",  sd, ',')



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
writedlm("se_scan77",  se, ',')



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
writedlm("he_scan77",  he, ',')


current_matrices = copy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)



for i in 0:0.005:1
    global prot_expression_cost
    current_Parameters["s_a"] = i
    current_Parameters["s_c"] = i
    current_matrices["selection_matrix"] = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))
    sup, dur, Y, S = timecourse_parameterscan_comprehensive(150, genotype_standard)
    add = [i sup dur Y S]
    prot_expression_cost = vcat(prot_expression_cost, add)

end
writedlm("prot_expression_cost_scan77",  prot_expression_cost, ',')



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
writedlm("sm_scan77",  sm, ',')


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
writedlm("ee_scan77",  ee, ',')



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
writedlm("eh_scan77",  eh, ',')




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
writedlm("es_scan77",  es, ',')


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

writedlm("he2_scan77",  he2, ',')

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

writedlm("shredding_activity_cost_scan77",  shredding_activity_cost, ',')


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

writedlm("m1_scan77",  m1, ',')

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

writedlm("m2_scan77",  m2, ',')

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

writedlm("er1_scan77",  er1, ',')



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
writedlm("er2_scan77",  er2, ',')


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
writedlm("er3_scan77",  er3, ',')



# ------------------
# also Rm
#
# apply_parameters_set(Parameters_set_sensitivity6_highRM)
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
# writedlm("Rm_scan76",  Rm, ',')
#
#
#
# # -------------------
# apply_parameters_set(Parameters_set_sensitivity7_highRM)
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
# writedlm("Rm_scan77",  Rm, ',')
