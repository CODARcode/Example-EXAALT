
from codar.cheetah import Campaign
from codar.cheetah import parameters as p

class Exaalt(Campaign):
	name = "Exaalt"
	codes = dict(exaalt="pt_producer_global",
				 stage_write="./stage_write/stage_write")

	supported_machines = ['titan_fob']

	project = "CSC143"
	queue = "batch"

	inputs = ["states_list.txt"]
	sweeps = [
		p.SweepGroup (nodes=80,
			post_processing = "",

			parameter_groups = 
			[
				p.Sweep([

					p.ParamRunner("exaalt", "nprocs", [1024]),
					p.ParamCmdLineArg("exaalt", "states_list_file", 1, ["states_list.txt"]),
					p.ParamCmdLineArg("exaalt", "no_of_states", 2, [1500000]),
					p.ParamCmdLineArg("exaalt", "bp_output_file", 3, ["output.bp"]),
					p.ParamCmdLineArg("exaalt", "transport_method", 4, ["FLEXPATH"]),
					p.ParamCmdLineArg("exaalt", "transport_variables", 5, [""]),
					p.ParamCmdLineArg("exaalt", "transport_options", 6, ["none"]),

					p.ParamRunner("stage_write", "nprocs", [256]),
					p.ParamCmdLineArg("stage_write", "input_bp_file", 1, ["output.bp"]),
					p.ParamCmdLineArg("stage_write", "output_bp_file", 2, ["staged.bp"]),
					p.ParamCmdLineArg("stage_write", "adios_read_method", 3, ["FLEXPATH"]),
					p.ParamCmdLineArg("stage_write", "read_method_params", 4, [""]),
					p.ParamCmdLineArg("stage_write", "adios_write_method", 5, ["POSIX"]),
					p.ParamCmdLineArg("stage_write", "write_method_params", 6, ["have_metadata_file=0"]),
					p.ParamCmdLineArg("stage_write", "variables_to_transform", 7, ["atom_id,atom_type,px,py,pz,imx,imy,imz,atom_vid,vx,vy,vz"]),
					p.ParamCmdLineArg("stage_write", "transform_params", 8, ["none"]),
				]),
			]),

		p.SweepGroup (nodes=80,
			post_processing = "",

			parameter_groups = 
			[
				p.Sweep([

					p.ParamRunner("exaalt", "nprocs", [1024]),
					p.ParamCmdLineArg("exaalt", "states_list_file", 1, ["states_list.txt"]),
					p.ParamCmdLineArg("exaalt", "no_of_states", 2, [1500000]),
					p.ParamCmdLineArg("exaalt", "bp_output_file", 3, ["output.bp"]),
					p.ParamCmdLineArg("exaalt", "transport_method", 4, ["FLEXPATH"]),
					p.ParamCmdLineArg("exaalt", "transport_variables", 5, [""]),
					p.ParamCmdLineArg("exaalt", "transport_options", 6, ["none"]),

					p.ParamRunner("stage_write", "nprocs", [256]),
					p.ParamCmdLineArg("stage_write", "input_bp_file", 1, ["output.bp"]),
					p.ParamCmdLineArg("stage_write", "output_bp_file", 2, ["staged.bp"]),
					p.ParamCmdLineArg("stage_write", "adios_read_method", 3, ["FLEXPATH"]),
					p.ParamCmdLineArg("stage_write", "read_method_params", 4, [""]),
					p.ParamCmdLineArg("stage_write", "adios_write_method", 5, ["POSIX"]),
					p.ParamCmdLineArg("stage_write", "write_method_params", 6, ["have_metadata_file=0"]),
					p.ParamCmdLineArg("stage_write", "variables_to_transform", 7, ["atom_id,atom_type,px,py,pz,imx,imy,imz,atom_vid,vx,vy,vz"]),
					p.ParamCmdLineArg("stage_write", "transform_params", 8, ["zlib:9"]),
				]),
			]),
		p.SweepGroup (nodes=80,
			post_processing = "",

			parameter_groups = 
			[
				p.Sweep([

					p.ParamRunner("exaalt", "nprocs", [1024]),
					p.ParamCmdLineArg("exaalt", "states_list_file", 1, ["states_list.txt"]),
					p.ParamCmdLineArg("exaalt", "no_of_states", 2, [1500000]),
					p.ParamCmdLineArg("exaalt", "bp_output_file", 3, ["output.bp"]),
					p.ParamCmdLineArg("exaalt", "transport_method", 4, ["FLEXPATH"]),
					p.ParamCmdLineArg("exaalt", "transport_variables", 5, [""]),
					p.ParamCmdLineArg("exaalt", "transport_options", 6, ["none"]),

					p.ParamRunner("stage_write", "nprocs", [256]),
					p.ParamCmdLineArg("stage_write", "input_bp_file", 1, ["output.bp"]),
					p.ParamCmdLineArg("stage_write", "output_bp_file", 2, ["staged.bp"]),
					p.ParamCmdLineArg("stage_write", "adios_read_method", 3, ["FLEXPATH"]),
					p.ParamCmdLineArg("stage_write", "read_method_params", 4, [""]),
					p.ParamCmdLineArg("stage_write", "adios_write_method", 5, ["POSIX"]),
					p.ParamCmdLineArg("stage_write", "write_method_params", 6, ["have_metadata_file=0"]),
					p.ParamCmdLineArg("stage_write", "variables_to_transform", 7, ["atom_id,atom_type,px,py,pz,imx,imy,imz,atom_vid,vx,vy,vz"]),
					p.ParamCmdLineArg("stage_write", "transform_params", 8, ["bzip2:9"]),
				]),
			]),
		p.SweepGroup (nodes=64,
			post_processing = "",

			parameter_groups = 
			[
				p.Sweep([

					p.ParamRunner("exaalt", "nprocs", [1024]),
					p.ParamCmdLineArg("exaalt", "states_list_file", 1, ["states_list.txt"]),
					p.ParamCmdLineArg("exaalt", "no_of_states", 2, [1500000]),
					p.ParamCmdLineArg("exaalt", "bp_output_file", 3, ["output.bp"]),
					p.ParamCmdLineArg("exaalt", "transport_method", 4, ["POSIX"]),
					p.ParamCmdLineArg("exaalt", "transport_variables", 5, ["have_metadata_file=0"]),
					p.ParamCmdLineArg("exaalt", "transport_options", 6, ["none"]),
					
				]),
			]),
	]
