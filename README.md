# 分子结构优化模型汇总
## 模型
### AutoGrow4
- [Paper: AutoGrow4: An open-source genetic algorithm for de novo drug design and lead optimization](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7165399/)
- [Code:](https://durrantlab.pitt.edu/autogrow4/)
- [论文解析](https://blog.csdn.net/weixin_42486623/article/details/131392519)
- [源码解析](https://blog.csdn.net/weixin_42486623/article/details/131409905)

script:
```bash
python RunAutogrow.py \
--filename_of_receptor tutorial/PARP/4r6eA_PARP1_prepared.pdb \
--center_x -70.76 \
--center_y 21.82 \
--center_z 28.33 \
--size_x 25.0 \
--size_y 16.0 \
--size_z 25.0 \
--source_compound_file source_compounds/naphthalene_smiles.smi \
--root_output_folder ./ \
--number_of_mutants_first_generation 50 \
--number_of_crossovers_first_generation 50 \
--number_of_mutants 50 \
--number_of_crossovers 50 \
--top_mols_to_seed_next_generation 50 \
--number_elitism_advance_from_previous_gen 50 \
--number_elitism_advance_from_previous_gen_first_generation 10 \
--diversity_mols_to_seed_first_generation 10  \
--diversity_seed_depreciation_per_gen 10 \
--num_generations 5 \
--mgltools_directory /DYZ/dyz1/download/mgltools_x86_64Linux2_1.5.7 \
--number_of_processors 1 \
--scoring_choice VINA \
--LipinskiLenientFilter \
--start_a_new_run \
--rxn_library click_chem_rxns \
--selector_choice Rank_Selector \
--dock_choice VinaDocking \
--max_variants_per_compound 5 \
--redock_elite_from_previous_gen False \
--generate_plot True \
--reduce_files_sizes True \
--use_docked_source_compounds True
```
