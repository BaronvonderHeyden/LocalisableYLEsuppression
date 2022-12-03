er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]


apply_parameters_set(Parameters_set_sensitivity2)

save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity2)

current_matrices = deepcopy(save_current_matrices)
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

writedlm("Simulation results\\er1_scan_005",  er1, ',')



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
writedlm("Simulation results\\er2_scan_005",  er2, ',')



#--------------------------------------
#--------------------------------------

er1 = Any["value" "sup" "dur" "Y" "S"]
er2 = Any["value" "sup" "dur" "Y" "S"]
er3 = Any["value" "sup" "dur" "Y" "S"]



apply_parameters_set(Parameters_set_sensitivity3)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity3)

current_matrices = deepcopy(save_current_matrices)
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

writedlm("Simulation results\\er1_scan_06",  er1, ',')



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
writedlm("Simulation results\\er2_scan_06",  er2, ',')


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
writedlm("Simulation results\\er3_scan",  er3, ',')

