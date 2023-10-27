
%%%% Autor: Philipp Amrein, University Freiburg, Medical Center, Radiology,
%%%% Medical Physics
%%%% February 2022

%This script generates a "x" gradient on a biplanar rectangular support
%strcuture

clc; clear all; close all; 

if ispc
cd('..\');
else
cd('../');
end

%% Run the algorithm


 coil_layouts.out=CoilGen(...
    'coil_mesh_file','bi_planer_rectangles_width_90.9mm_distance_32.85mm.stl', ...
    'field_shape_function','y',...
    'secondary_target_mesh_file','none', ...
    'secondary_target_weight',0, ...
    'target_region_radius',0.0085,...  % in meter
    'conductor_thickness', 0.00058,...
    'use_only_target_mesh_verts',false, ...
    'sf_source_file','none', ...
    'levels',10, ... % the number of potential steps that determines the later number of windings (Stream function discretization)
    'pot_offset_factor',0.5, ... % a potential offset value for the minimal and maximal contour potential ; must be between 0 and 1
    'surface_is_cylinder_flag',false, ...
    'interconnection_cut_width',0.0058, ... % the width for the interconnections are interconnected; in meter
    'normal_shift_length',0.001, ... % the length for which overlapping return paths will be shifted along the surface normals; in meter
    'iteration_num_mesh_refinement',1, ... % the number of refinements for the mesh;
    'set_roi_into_mesh_center',true, ...
    'conductor_cross_section_height',0.00003, ...
    'conductor_cross_section_width',0.00003, ...
    'force_cut_selection',{'high'},...
    'level_set_method','primary',... %Specify one of the three ways the level sets are calculated: "primary","combined", or "independent"
    'interconnection_method','regular',...
    'skip_postprocessing',false,...
    'skip_inductance_calculation',true,...
    'tikonov_reg_factor',1000); %Tikonov regularization factor for the SF optimization


%% Plot results
close all;

coil_name='Coil';

if ispc
addpath(strcat(pwd,'\','plotting'));
else
addpath(strcat(pwd,'/','plotting'));
end
%Chose a even leveled solution for plotting
solutions_to_plot=find(arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)));
single_ind_to_plot= find_even_leveled_solution(coil_layouts);
plot_error_different_solutions(coil_layouts,single_ind_to_plot,coil_name);
plot_2D_contours_with_sf(coil_layouts,single_ind_to_plot,coil_name);
plot_groups_and_interconnections(coil_layouts,single_ind_to_plot,coil_name);
plot_coil_parameters(coil_layouts,coil_name);
plot_coil_track_with_resulting_bfield(coil_layouts,single_ind_to_plot,coil_name);
%plot_various_error_metrics(coil_layouts,single_ind_to_plot,coil_name);
plot_resulting_gradient(coil_layouts,single_ind_to_plot,coil_name);
rmpath('plotting');

