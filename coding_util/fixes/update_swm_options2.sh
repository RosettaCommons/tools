# if you are on a mac, change "-i" to "-i .old" below.
for i in "$@"; do 

sed -i '~' 's/verbose_scores_/options_->verbose_scores()/g;s/force_centroid_interaction_/options_->force_centroid_interaction()/g;s/use_phenix_geo_/options_->use_phenix_geo()/g;s/skip_deletions_/options_->skip_deletions()/g;s/erraser_/options_->erraser()/g;s/allow_internal_moves_/options_->allow_internal_moves()/g;s/num_random_samples_/options_->num_random_samples()/g;s/cycles_/options_->cycles()/g;s/add_delete_frequency_/options_->add_delete_frequency()/g;s/minimize_single_res_frequency_/options_->minimize_single_res_frequency()/g;s/minimize_single_res_/options_->minimize_single_res()/g;s/minimizer_allow_variable_bond_geometry_/options_->minimizer_allow_variable_bond_geometry()/g;s/minimizer_vary_bond_geometry_frequency_/options_->minimizer_vary_bond_geometry_frequency()/g;s/switch_focus_frequency_/options_->switch_focus_frequency()/g;s/just_min_after_mutation_frequency_/options_->just_min_after_mutation_frequency()/g;s/temperature_/options_->temperature()/g;s/sample_res_/options_->sample_res()/g;s/bulge_res_/options_->bulge_res()/g;s/allow_skip_bulge_/options_->allow_skip_bulge()/g;s/virtual_sugar_keep_base_fixed_/options_->virtual_sugar_keep_base_fixed()/g;s/constraint_x0_/options_->constraint_x0()/g;s/constraint_tol_/options_->constraint_tol()/g;s/extra_minimize_res_/options_->extra_minimize_res()/g;s/syn_chi_res_list_/options_->syn_chi_res_list()/g' $i
done




