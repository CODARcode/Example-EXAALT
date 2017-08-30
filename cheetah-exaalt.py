
from codar.cheetah import Campaign
from codar.cheetah import parameters as p

class Exaalt(Campaign):
	name = "Exaalt"
	codes = dict(exaalt="pt_producer_global",
				 stage_write="../heat-transfer-codar/stage_write/stage_write")

	supported_machines = ['titan_fob']

	project = "CSC143"
	queue = "debug"

	inputs = ["states_list.txt"]
	sweeps = [
		p.SweepGroup (nodes=20,
			post_processing = "",

			parameter_groups = 
			[
				p.Sweep([

					p.ParamRunner("exaalt", "nprocs", [8,16]),
					p.ParamCmdLineArg("exaalt", "states_list_file", 1, ["states_list.txt"]),
					p.ParamCmdLineArg("exaalt", "no_of_states", 2, [16]),
					p.ParamCmdLineArg("exaalt", "no_of_randoms", 3, [16]),
					p.ParamCmdLineArg("exaalt", "bp_output_file", 4, ["output.bp"]),
					p.ParamCmdLineArg("exaalt", "transport_method", 5, ["FLEXPATH"]),
					p.ParamCmdLineArg("exaalt", "transport_options", 6, [""]),

					p.ParamRunner("stage_write", "nprocs", [2]),
					p.ParamCmdLineArg("stage_write", "input_bp_file", 1, ["output.bp"]),
					p.ParamCmdLineArg("stage_write", "output_bp_file", 2, ["staged.bp"]),
					p.ParamCmdLineArg("stage_write", "adios_read_method", 3, ["FLEXPATH"]),
					p.ParamCmdLineArg("stage_write", "read_method_params", 4, [""]),
					p.ParamCmdLineArg("stage_write", "adios_write_method", 5, ["MPI"]),
					p.ParamCmdLineArg("stage_write", "write_method_params", 6, [""]),
					p.ParamCmdLineArg("stage_write", "variables_to_transform", 7, ["px,py,pz"]),
					p.ParamCmdLineArg("stage_write", "transform_params", 8, ["none","zlib", "bzip2"]),
					p.ParamCmdLineArg("stage_write", "decomposition", 9, [2]),
				]),
			]),
	]
