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


    end



    if length(network_ids) > 1
        # Initialize head variables for the final time index.
        variable_head(wm; nw = network_ids[end])

        # Constraints on network expansion variables.
        # for n_2 in network_ids[2:end-1]
        #     # Constrain short pipe selection variables based on the initial time index.
        #     for i in ids(wm, :ne_short_pipe; nw = n_2)
        #         constraint_ne_short_pipe_selection(wm, i, n_1, n_2)
        #     end
        # end
    end

    objective_min_ec_ls(wm)
end
