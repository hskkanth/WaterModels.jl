# function solve_ed(network, model_constructor, optimizer; kwargs...)
#     return solve_model(network, model_constructor, optimizer, build_ed; kwargs...)
# end


function solve_mn_ec_ls(file, model_constructor, optimizer; kwargs...)
    return solve_model(file, model_constructor, optimizer, build_mn_ec_ls; multinetwork=true, kwargs...)
end


function objective_min_ec_ls(wm::AbstractWaterModel)::JuMP.AffExpr

    println("Check scaling")

    # Get all network IDs in the multinetwork.
    first_network_id = sort(collect(nw_ids(wm)))[1]

    # Initialize the objective expression to zero.
    objective = JuMP.AffExpr(0.0)
  
    n = first_network_id
    # Get the set of network expansion short pipes at time index `n`.
    for (a, ne_short_pipe) in ref(wm, n, :ne_short_pipe)
        # Add the cost of network expansion component `a` at time period `n`.
        term = ne_short_pipe["construction_cost"] * var(wm, n, :z_ne_short_pipe, a)
        JuMP.add_to_expression!(objective, term)
    end

    # Get the set of design pipes at time index `n`.
    for (a, des_pipe) in ref(wm, n, :des_pipe)
        # Add the cost of network expansion component `a` at time period `n`.
        term = des_pipe["construction_cost"] * var(wm, n, :z_des_pipe, a)
        JuMP.add_to_expression!(objective, term)
    end

    # Get the set of network expansion pumps at time index `n`.
    for (a, ne_pump) in ref(wm, n, :ne_pump)
        # Add the cost of network expansion component `a` at time period `n`.
        term = ne_pump["construction_cost"] * var(wm, n, :x_ne_pump, a)
        JuMP.add_to_expression!(objective, term)
    end

    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    # Find the network IDs to define the objective over.
    if length(network_ids) > 1
        network_ids_flow = network_ids[1:end-1]
    else
        network_ids_flow = network_ids
    end

    for n in network_ids_flow
        # Get the set of dispatchable demands at time index `n`.
        dispatchable_demands = ref(wm, n, :dispatchable_demand)

        for (i, demand) in filter(x -> x.second["flow_min"] >= 0.0, dispatchable_demands)
            # Add the volume delivered at demand `i` and time period `n`.
            
            base_time = wm.ref[:it][wm_it_sym][:base_time]
            volume_conversion = get(demand, "priority", 1.0) * ref(wm, n, :time_step) * base_time
            shedding_cost = 1.0;
            base_flow = wm.ref[:it][wm_it_sym][:base_flow]
            JuMP.add_to_expression!(objective, shedding_cost * volume_conversion * base_flow * (demand["flow_max"]- var(wm, n, :q_demand, i)))
        end
    end

    # Minimize the total cost of network expansion.
    return JuMP.@objective(wm.model, JuMP.MIN_SENSE, objective)
end

function build_mn_ec_ls(wm::AbstractWaterModel)
    # Get all network IDs in the multinetwork.
    network_ids = sort(collect(nw_ids(wm)))

    if length(network_ids) > 1
        network_ids_inner = network_ids[1:end-1]
    else
        network_ids_inner = network_ids
    end

    # Start with the first network, representing the initial time step.
    n_1 = network_ids[1]
    # Variables
        # Binary variables
        variable_ne_short_pipe_indicator(wm; nw=n_1)
        variable_des_pipe_indicator(wm; nw=n_1)
        variable_ne_pump_build(wm; nw=n_1)
    # Components involved: reservoirs, pipes, demands
    for n in network_ids_inner

        # Main continuous variables
        variable_head(wm; nw = n)
        variable_flow(wm; nw = n)

        # Node Attachment variables
        variable_demand_flow(wm; nw = n)
        variable_reservoir_flow(wm; nw = n)
        variable_tank_flow(wm; nw=n)

        #Pump
        variable_pump_head_gain(wm; nw=n)
        variable_pump_indicator(wm; nw=n)
        # variable_pump_power(wm; nw=n)
        #Ne Pump
        variable_ne_pump_head_gain(wm; nw=n)
        variable_ne_pump_indicator(wm; nw=n)
        # variable_ne_pump_power(wm;nw=n)

    #Constraints
        # Flow conservation at all nodes.
        for i in ids(wm, :node; nw = n)
            constraint_flow_conservation(wm, i; nw = n)
            constraint_node_directionality(wm, i; nw = n)
        end
        # Constraints on pipe flows, heads, and physics.
        for a in ids(wm, :pipe; nw=n)
            constraint_pipe_flow(wm, a; nw=n)
            constraint_pipe_head(wm, a; nw=n)
            constraint_pipe_head_loss(wm, a; nw=n)
        end

        # # Constraints on design pipe flows, heads, and physics.
        for a in ids(wm, :des_pipe; nw=n)
            constraint_on_off_des_pipe_head(wm, a; nw=n)
            constraint_on_off_des_pipe_head_loss(wm, a; nw=n)
            constraint_on_off_des_pipe_flow(wm, a; nw=n)
        end

        # Constraints on short pipe flows and heads.
        for a in ids(wm, :short_pipe; nw=n)
            constraint_short_pipe_head(wm, a; nw=n)
            constraint_short_pipe_flow(wm, a; nw=n)
        end

        # Constraints on expansion short pipe flows and heads.
        for a in ids(wm, :ne_short_pipe; nw=n)
            constraint_short_pipe_head_ne(wm, a; nw=n)
            constraint_short_pipe_flow_ne(wm, a; nw=n)
        end

        # Constraints on pump flows, heads, and physics.
        for a in ids(wm, :pump; nw=n)
            constraint_on_off_pump_head(wm, a; nw=n)
            constraint_on_off_pump_head_gain(wm, a; nw=n)
            constraint_on_off_pump_flow(wm, a; nw=n)
            # constraint_on_off_pump_power(wm, a; nw=n)
        end

        # Constraints on expansion pump flows, heads, and physics.
        for a in ids(wm, :ne_pump; nw=n)
            constraint_on_off_pump_head_ne(wm, a; nw=n)
            constraint_on_off_pump_head_gain_ne(wm, a; nw=n)
            constraint_on_off_pump_flow_ne(wm, a; nw=n)
            # constraint_on_off_pump_power_ne(wm, a; nw=n)
            constraint_on_off_pump_build_ne(wm, a; nw=n)
        end
    end


    multinetwork_tank_constraints(wm, network_ids)
    
    
    ###### End of Tank constraints for multi-network formulation ######

    objective_min_ec_ls(wm)
end


function multinetwork_tank_constraints(wm::AbstractWaterModel, network_ids::Vector{Int64})
     ###### Tank constraints for multi-network formulation ######
    # Set initial conditions of tanks.
    n_1 = network_ids[1]
    n_f = network_ids[end]

    tank_volume_reset_time_points = Set(get(wm.ref[:it][wm_it_sym], :tank_volume_reset_time_points, Int[]))
    tank_volume_recovery_time_points = Set(get(wm.ref[:it][wm_it_sym], :tank_volume_recovery_time_points, Int[]))

    push!(tank_volume_reset_time_points, n_1)
    push!(tank_volume_recovery_time_points, n_f)


    if length(network_ids) > 1
        
        # Initialize head variables for the final time index.
        variable_head(wm; nw = network_ids[end])

        n_prev = n_1
        # Constraints on tank volumes.
        for n in network_ids
            for i in ids(wm, :tank; nw = n)
                if n in tank_volume_reset_time_points
                    @info "If pass for reset; n = $n, n_prev = $n_prev"
                    constraint_tank_volume(wm, i; nw = n) #includes n_1
                else
                    constraint_tank_volume(wm, i, n_prev, n)
                end
            end

            # Update the first network used for integration.
            n_prev = n
        end

        # Ensure tanks recover their initial volume.
        for n_tank in tank_volume_recovery_time_points
            for i in ids(wm, n_tank, :tank)
                constraint_tank_volume_recovery(wm, i, n_1, n_tank)
            end
        end

    else #single time period case;
        for i in ids(wm, :tank; nw = n_1)
            constraint_tank_volume(wm, i; nw = n_1)
        end
    end
end


 