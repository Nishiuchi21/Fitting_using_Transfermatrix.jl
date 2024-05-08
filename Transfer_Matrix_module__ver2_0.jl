module Transfer_Matrix_module

    using GLMakie
    using TransferMatrix
    using DataFrames
    using CSV
    using Optim
    using LineSearches
    using RefractiveIndex



    function wavenumber2wavelength(ν)
        return 1.0 ./ ν
    end

    function wavelength2wavenumber(λ)
        return 1.0 ./ λ
    end

    function wavelength_range(ν_top, ν_bottom, resolution)
        λs = range(wavenumber2wavelength(ν_top * 1e2), wavenumber2wavelength(ν_bottom * 1e2), length = resolution)
        println("Wavelength range: $(round(collect(λs)[1] * 1e6, digits = 2)) - $(round(collect(λs)[end] * 1e6, digits = 2)) (μm)")
        return λs
    end

    function initialize_figure(resolution::Tuple, fontsize::Int; data_inspector = true)
        f = Figure(resolution = resolution,fontsize = fontsize)
        display(f)
        if data_inspector == true
            DataInspector(f)
        end

        return f
    end

    function initialize_axis(gridpos::GridPosition ,title::String, titlesize::Int, xlabel::String, ylabel::String, label_size::Int, xticks, yticks; second_axis = false, second_axis_yticks = LinearTicks(5), second_axis_ylabel = "")
        ax_list = []
        ax = Axis(gridpos, title = title, titlesize = titlesize, xlabel = xlabel, xlabelsize = label_size, xlabelfont = :bold, ylabel = ylabel, ylabelsize = label_size, ylabelfont = :bold, xticks = xticks, yticks = yticks)
        push!(ax_list, ax)
        if second_axis == true
            ax_second = Axis(gridpos, yaxisposition = :right, ylabel = second_axis_ylabel, ylabelsize = label_size, ylabelfont = :bold, yticks = second_axis_yticks)
            hidexdecorations!(ax_second)
            hidedecorations!(ax_second, label = false, ticklabels = false, ticks = false, minorticks = false)
            push!(ax_list, ax_second)
        end

        return ax_list
    end

    function total_height(s)
        height = 0
        for i in 1:length(s.layers)
            height += s.layers[i].thickness
        end

        return height
    end

    function height_list_creator(s)
        height_list = []
        for i in 1:length(s.layers)
            height = - s.layers[1].thickness
            for j in 1:(i - 1)
                height += s.layers[j].thickness
            end

            push!(height_list, height + (s.layers[i].thickness))
        end

        return height_list
    end

    function plot_band(s, color_dict::Dict, ax)
        hight_list = hight_list_creator(s) * 1e6
        pushfirst!(hight_list, - s.layers[1].thickness * 1e6)
        println(ax)

        for i in 1:length(s.layers)
            color = color_dict[s.layers[i].material]
            band!(range(3.5, 4.5), hight_list[i], hight_list[i + 1], color = color) 
            text!(ax, s.layers[i].material * " : " * "$(round(s.layers[i].thickness * 1e6, digits = 4))" * "(μm)", position = (4.6, hight_list[i]))
        end
    end

    function initialize_legend(ax::Axis, All_lines, legend_names::Array; position = (0.01, 0.99), patch_size = (30,50), label_size = 30)
        legend = axislegend(ax, All_lines, legend_names ,position = position, patchsize = patch_size,labelsize = label_size)
        return legend
    end

    function choose_files(l::Array, index_array::Array)
        l_plot = []
        for i in index_array
            push!(l_plot, l[i])
        end

        return l_plot
    end

    function read_folder(root::String)
        l = readdir(root; join = true)
        # println("\nFiles in the folder: ")
        # for (i, j) in enumerate(l)
        #     println("No.$(i): $(j)")
        # end

        return l
    end

    function peaks_file_reader(file)
        data = DataFrame(CSV.File(abspath(file)))
        rename!(data,["Column1"])

        return data.Column1
    end


    function actuatable_DBR_cavity(cav_length, λs; refractive_index_correction_zns = -0.07709, refractive_index_correction_ge = 0.03968, amp = 0.8309, material_refractive_index = 1.0)
        
        l_cav = cav_length * 1e-6
        material_r = material_refractive_index    
        material = Layer("Material",l_cav,collect(λs),fill((material_r),length(λs)),zeros(Float64,length(λs)))


        zns_autodata = RefractiveMaterial("main","ZnS","Amotchkina") 
        d_zns = 515.12 * 1e-9
        
        ge_autodata = RefractiveMaterial("main","Ge","Icenogle")
        d_ge = 290.66 * 1e-9
        
        caf2_autodata = RefractiveMaterial("main","CaF2","Malitson")
        d_caf2 = 3e-3
        
        zns = Layer(zns_autodata[1].name, d_zns, λs, zns_autodata[1].(λs .* 1e6), zeros(length(λs)))
        ge = Layer(ge_autodata.name, d_ge, λs, ge_autodata.(λs .* 1e6), zeros(length(λs)))
        caf2 = Layer(caf2_autodata.name, d_caf2, λs, caf2_autodata.(λs .* 1e6), zeros(length(λs)))
        
        zns.n .+= refractive_index_correction_zns
        ge.n .+= refractive_index_correction_ge


        layers_1 = Layer[]
        n_repeat = 4
        for i in 1:n_repeat
            push!(layers_1,zns)
            push!(layers_1,ge)
        end
        push!(layers_1,caf2)

        layers_2 = reverse(layers_1)
        layers = vcat(layers_2,[material],layers_1)


        s = Structure(layers, collect(λs), [0.0])
        Tp,Ts,Rp,Rs = calculate_tr(s)

        return Tp .* amp
    end



    function single_DBR_mirror(λs; refractive_index_correction_zns = -0.07709, refractive_index_correction_ge = 0.03968, amp = sqrt(0.8309))


        zns_autodata = RefractiveMaterial("main","ZnS","Amotchkina") 
        d_zns = 515.12 * 1e-9
        
        ge_autodata = RefractiveMaterial("main","Ge","Icenogle")
        d_ge = 290.66 * 1e-9
        
        caf2_autodata = RefractiveMaterial("main","CaF2","Malitson")
        d_caf2 = 3e-3
        
        zns = Layer(zns_autodata[1].name, d_zns, λs, zns_autodata[1].(λs .* 1e6), zeros(length(λs)))
        ge = Layer(ge_autodata.name, d_ge, λs, ge_autodata.(λs .* 1e6), zeros(length(λs)))
        caf2 = Layer(caf2_autodata.name, d_caf2, λs, caf2_autodata.(λs .* 1e6), zeros(length(λs)))
        
        zns.n .+= refractive_index_correction_zns
        ge.n .+= refractive_index_correction_ge


        layers_1 = Layer[]
        n_repeat = 4
        for i in 1:n_repeat
            push!(layers_1,zns)
            push!(layers_1,ge)
        end
        push!(layers_1,caf2)

        layers_2 = reverse(layers_1)
        layers = layers_2


        s = Structure(layers, collect(λs), [0.0])
        Tp,Ts,Rp,Rs = calculate_tr(s)

        return Tp .* amp
    end



    function peak_finder(df_res ; threshold = 30)
        peak_position_list = []
        for i in 2:(length(df_res.wavenumber) - 1)
            if df_res.T[i] > threshold
                if df_res.T[i] >= df_res.T[i + 1] && df_res.T[i] >= df_res.T[i - 1]
                    peak_position = df_res.wavenumber[i]
                    push!(peak_position_list,peak_position)
                end    
            end
        end

        return peak_position_list
    end
    
    
    
    function error(p, x, y; material_r = 1.0)
        error = 0

        df = DataFrame(wavenumber = x, T = reverse!(actuatable_DBR_cavity(p[1], reverse!(wavenumber2wavelength.(x * 1e2)), material_refractive_index = material_r) .* 100))

        for i in 1:length(x)
            error += (y[i] - df.T[i])^2
        end

        println("error: ", error)
        println("cavity length:", p[1])
        return error
    end



    function error_resonance_cavity_length(p, calculation_range, Abs_peak, material_r)
        error = 0

        df = DataFrame(wavenumber = calculation_range, T = reverse!(Transfer_Matrix_module.actuatable_DBR_cavity(p[1], reverse!(Transfer_Matrix_module.wavenumber2wavelength.(calculation_range .* 1e2)), material_refractive_index = material_r) .* 100))
        # println(df.T[1:10])
        peaks = peak_finder(df, threshold = 60)
        # println(peaks)

        error = (peaks[1] - Abs_peak)^2

        println("error: ", error)
        println("cavity length:", p[1])
        return error
    end



    function error_mt(p, x, y) # multi-threaded. However, it doesn't make the calculation faster due to the low number of calculations.
        df = DataFrame(wavenumber = x, T = reverse!(actuatable_DBR_cavity(p[1], reverse!(wavenumber2wavelength.(x * 1e2)), refractive_index_correction_zns = p[2]) .* 100))
        
        
        error = zeros(Float64, Threads.nthreads())
        Threads.@threads for i in 1:length(x)
            error[Threads.threadid()] += (y[i] - df.T[i])^2
        end
        sum_error = sum(error)

        println("error: ", sum_error)
        println("cavity length:", p[1])
        return sum_error
    end


    function graph_data(file; option=false,range = [1000,2000])
        data = DataFrame(CSV.File(abspath(file)))
        rename!(data,["Column1","Column2"])

        if data.Column2[1] == "INFRARED SPECTRUM"
            data = DataFrame(CSV.File(abspath(file),skipto = 20,footerskip = 37))
            println("\n $(file) \n !!!This data is from the FTIR spectrometer!!!")
        end

        rename!(data,["wavenumber","transmittance"])
        
        min,max = range
        if option == true
            data = data[(data.wavenumber .< max) .& (data.wavenumber .> min),:]
        end
    
        return data.wavenumber,data.transmittance
    end

    function find_cavity_length(l::Vector, index::Int,  fitting_parameters::Vector, material_r::Float64; range = [1800, 2400], div_factor = 1, rt = 0.8, min_bounds = [0.0], max_bounds = [6.0])
        data = l[index]
        x,y = graph_data(data, option = true, range = range)
        x = [x[div_factor * i] for i in 1:Int(length(x) ÷ div_factor)]
        y = [y[div_factor * i] for i in 1:Int(length(y) ÷ div_factor)]

        x_whole,y_whole = graph_data(data)

        # result = optimize(b -> error(b, x, y), min_bounds, max_bounds, fitting_parameters, SAMIN(rt = rt, neps = 20))
        result = optimize(b -> error(b, x, y, material_r = material_r), fitting_parameters, NelderMead(parameters = Optim.FixedParameters(α=0.5, β=1., γ=0.8, δ=0.8))) # Probably, α and β are for rough searches, amd γ and δ are for fine searches.
        p = Optim.minimizer(result)

        return p, x_whole
    end
    


    function resonance_cavity_length_calculator(Abs_peak::Float64, calculation_range::StepRangeLen,  first_guess::Vector, material_r::Float64; rt = 0.1, min_bounds = [0.0], max_bounds = [6.0])

        result = optimize(b -> error_resonance_cavity_length(b, calculation_range, Abs_peak, material_r), min_bounds, max_bounds, first_guess, SAMIN(rt = rt, neps = 20))
        # result = optimize(b -> error_resonance_cavity_length(b, calculation_range, Abs_peak, material_r), first_guess, NelderMead(parameters = Optim.FixedParameters(α=0.5, β=1., γ=0.8, δ=0.8))) # Probably, α and β are for rough searches, amd γ and δ are for fine searches.
        p = Optim.minimizer(result)

        return p[1]
    end



    function plot(l_plot::Array , ax::Axis, linewidth::Int; linestyle = :solid, scatter = false)
        lines_l = []
        for (i,j) in enumerate(l_plot)
            println("\nPlotting: ",j)
            x,y = graph_data(j)
            # lines!(ax,x,y)
            if scatter == true
                m = scatter!(ax,x,y,linewidth = linewidth, color = Cycled(i), linestyle = linestyle)
            else
                m = lines!(ax,x,y,linewidth = linewidth, color = Cycled(i), linestyle = linestyle)
            end
            
            push!(lines_l,m)
        end

        println(lines_l)
        return lines_l
    end


end
