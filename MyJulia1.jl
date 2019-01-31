import Pkg
if ~haskey(Pkg.installed(), "PowerModels")
	Pkg.add("PowerModels");
end
if ~haskey(Pkg.installed(), "Ipopt")
	Pkg.add("Ipopt");
end
using PowerModels, Ipopt
using DelimitedFiles
include("result_process_PowerModel.jl");
include("data_process.jl");

function MyJulia1(conFile, inlFile, rawFile, ropFile, timeLimitSeconds, scoringMethod, networkModel)
	solverI = IpoptSolver();

	paths="Solution1.m";
	network = PowerModels.parse_file(paths);
	data_process(network,1);
	result = run_pf(network, ACPPowerModel, solverI);
	(Pg, Qg, V, Theta) = result_process_PowerModel_PF(result);

	function create_soln_bus_gen_matrices(network)
		nbus = length(network["bus"]);
		bus_sec = Array{Any}(undef, nbus, 4); 
		for i = 1:nbus
		   bus_sec[i, 1] = network["bus"][string(i)]["bus_i"];
		   bus_sec[i, 2] = network["bus"][string(i)]["vm"];
		   bus_sec[i, 3] = network["bus"][string(i)]["va"];
		   bus_sec[i, 4] = Float64(0);
		end
		for i = 1:length(network["shunt"])
		   # TODO: Convert the bs into the correct units - possibly need to multiply by 100 
		   bus_sec[network["shunt"][string(i)]["shunt_bus"], 4] = network["shunt"][string(i)]["bs"];
		end

		# TODO: Convert the p and q into the correct units - possibly need to multiply by 100 
		ngen = length(network["gen"]);
		gen_sec = Array{Any}(undef, ngen, 4)
		for i = 1:ngen
		   gen_sec[i, 1] = network["gen"][string(i)]["gen_bus"];
		   gen_sec[i, 2] = "'"*string(network["gen"][string(i)]["index"])*"'";
		   gen_sec[i, 3] = network["gen"][string(i)]["pg"];
		   gen_sec[i, 4] = network["gen"][string(i)]["qg"];
		end
		return (bus_sec, gen_sec)
	end

	function write_bus_gen_sec(fhandle, basecase)
	   bussec = basecase["bus_section"]
	   write(fhandle, "--bus section\n");
	   write(fhandle, "i, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n");
	   writedlm(fhandle, bussec, ",");
	   gensec = basecase["generator_section"]
	   ngen, _ = size(gensec)
	   for i = 1:ngen
	      gensec[i, 2] = "'"*string(gensec[i, 2])*"'"
	   end
	   write(fhandle, "--generator section\n");
	   write(fhandle, "i, id, p(MW), q(MVAR)\n");
	   writedlm(fhandle, gensec, ",");
	end

	function write_sol1_file(basecase, filename)
		open(filename, "w") do f
			write_bus_gen_sec(f, basecase);           
		end
	end

	#function write_sol2_file(bussec_array, gensec_array, con_label_array, delta_array, filename)
	function write_sol2_file(contingency, filename)
		open(filename, "w") do f
			for i =1:length(contingency)
				write(f, "--contingency\n");
				write(f, "label\n");
				write(f, "'"*contingency[string(i)]["label"]*"'"*"\n");
				write_bus_gen_sec(f, contingency[string(i)]);
				write(f, "--delta section\n");
				write(f, "delta(MW)\n");
				write(f, string(contingency[string(i)]["delta"]*"\n"));
			end
		end
	end

	(bussec, gensec) = create_soln_bus_gen_matrices(network);
	write_sol1_file(bussec, gensec, "Solution1.txt");

	paths="Solution2.m";
	network = PowerModels.parse_file(paths);
	data_process(network,1);
	result = run_pf(network, ACPPowerModel, solverI);
	(Pg, Qg, V, Theta) = result_process_PowerModel_PF(result);

	(bussec, gensec) = create_soln_bus_gen_matrices(network);
	bussec_array = Any[];
	push!(bussec_array, bussec);
	gensec_array = Any[];
	push!(gensec_array, gensec);
	con_label_array = ["LINE-6-12-BL"];
	delta_array = [0];
	write_sol2_file(bussec_array, gensec_array, con_label_array, delta_array, "Solution2.txt");	
end
