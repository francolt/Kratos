# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv
import numpy as np
np.random.seed(2021)

# Read parameters
with open("stochastic_optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"])
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
original_directory = os.getcwd()
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"

# Testing by
# 1) using the "json_output_process" & "json_check_process" within the structural analysis
# 2) additionally checking some process output
os.chdir(output_directory)

with open(optimization_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_iteration = float(last_line[0].strip())
    resulting_abs_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_iteration, 3)
    TestCase().assertAlmostEqual(resulting_abs_improvement, -12.5586, 3)

os.chdir(original_directory)


kratos_utilities.DeleteFileIfExisting("current_design.mdpa")
kratos_utilities.DeleteFileIfExisting("current_design.time")

# =======================================================================================================