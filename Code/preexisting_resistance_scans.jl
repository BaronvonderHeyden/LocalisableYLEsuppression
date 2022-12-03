


shredding_a = Any["value" "max_sup" "t_max_sup" "t_1e_2" "t_1e_3" "t_1e_4" "t_1e_5" "t_1e_6" "t_1e_7" "t_1e_8" "t_1e_9" "t_1e_10" "dur" "Y" "S"]
homing_a = Any["value" "max_sup" "t_max_sup" "t_1e_2" "t_1e_3" "t_1e_4" "t_1e_5" "t_1e_6" "t_1e_7" "t_1e_8" "t_1e_9" "t_1e_10" "dur" "Y" "S"]
shredding_b =Any["value" "max_sup" "t_max_sup" "t_1e_2" "t_1e_3" "t_1e_4" "t_1e_5" "t_1e_6" "t_1e_7" "t_1e_8" "t_1e_9" "t_1e_10" "dur" "Y" "S"]
homing_b = Any["value" "max_sup" "t_max_sup" "t_1e_2" "t_1e_3" "t_1e_4" "t_1e_5" "t_1e_6" "t_1e_7" "t_1e_8" "t_1e_9" "t_1e_10" "dur" "Y" "S"]


apply_parameters_set(Parameters_set_sensitivity2)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity2)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global shredding_a

    release_pop = generate_population(r2=i)
    release_pop[return_i(["AB" "ef" "CD" "cd"], genotypes_detailed)] = BigFloat("0.001")

    sup,t_max,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,dur,Y,S = timecourse_parameterscan_preexisting_resistance(1000, release_pop)
    add = [i sup t_max e_2 e_3 e_4 e_5 e_6 e_7 e_8 e_9 e_10 dur Y S]
    shredding_a = vcat(shredding_a, add)

end

writedlm("Simulation results\\shredding_005_1000_generations",  shredding_a, ',')


for i in 0:0.01:1
    global homing_a

    release_pop = generate_population(r3=i)
    release_pop[return_i(["AB" "ef" "CD" "cd"], genotypes_detailed)] = BigFloat("0.001")


    sup,t_max,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,dur,Y,S = timecourse_parameterscan_preexisting_resistance(1000, release_pop)
    add = [i sup t_max e_2 e_3 e_4 e_5 e_6 e_7 e_8 e_9 e_10 dur Y S]
    homing_a = vcat(homing_a, add)

end

writedlm("Simulation results\\homing_005_1000_generations",  homing_a, ',')


apply_parameters_set(Parameters_set_sensitivity3)
save_current_matrices = deepcopy(current_matrices)
save_current_Parameters = deepcopy(Parameters_set_sensitivity3)

current_matrices = deepcopy(save_current_matrices)
current_Parameters = deepcopy(save_current_Parameters)


for i in 0:0.01:1
    global shredding_b

    release_pop = generate_population(r2=i)
    release_pop[return_i(["AB" "ef" "CD" "cd"], genotypes_detailed)] = BigFloat("0.001")

    sup,t_max,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,dur,Y,S = timecourse_parameterscan_preexisting_resistance(1000, release_pop)
    add = [i sup t_max e_2 e_3 e_4 e_5 e_6 e_7 e_8 e_9 e_10 dur Y S]
    shredding_b = vcat(shredding_b, add)

end

writedlm("Simulation results\\shredding_060_1000_generations",  shredding_b, ',')


for i in 0:0.01:1
    global homing_b
    release_pop = generate_population(r3=i)
    release_pop[return_i(["AB" "ef" "CD" "cd"], genotypes_detailed)] = BigFloat("0.001")

    sup,t_max,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9,e_10,dur,Y,S = timecourse_parameterscan_preexisting_resistance(1000, release_pop)
    add = [i sup t_max e_2 e_3 e_4 e_5 e_6 e_7 e_8 e_9 e_10 dur Y S]
    homing_b = vcat(homing_b, add)

end

writedlm("Simulation results\\homing_060_1000_generations",  homing_b, ',')

