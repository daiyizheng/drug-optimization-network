"""
This function should contain all the info for executing the mutation functions
"""

import __future__

import random
import copy


import autogrow.operators.mutation.smiles_click_chem.smiles_click_chem as SmileClickClass


#######################################
# Functions for creating molecular models
##########################################
def make_mutants(vars, generation_num, number_of_processors,
                 num_mutants_to_make, ligands_list,
                 new_mutation_smiles_list,
                 rxn_library_variables):
    """
    Make mutant compounds in a list to be returned
    在要返回的列表中制作突变化合物
    This runs SmileClick and returns a list of new molecules
    这将运行 SmileClick 并返回新分子列表
    Inputs:
    :param dict vars: a dictionary of all user variables
    :param int generation_num: generation number
    :param int number_of_processors: number of processors as specified by the
        user
    :param int num_mutants_to_make: number of mutants to return
    :param list ligands_list: list of ligand/name pairs which are the order in
        which to be sampled :param list new_mutation_smiles_list: is the list of
        mutants made for the current generation being populated but in a previous
        iteration of the loop in Operations
    :param list rxn_library_variables: a list of user variables which define
        the rxn_library, rxn_library_file, and function_group_library. ie.
        rxn_library_variables = [vars['rxn_library'], vars['rxn_library_file'],
        vars['function_group_library']]

    Returns:
    :returns: list new_ligands_list: ligand/name pairs OR returns None if
        sufficient number was not generated. None: bol if mutations failed
    """

    if len(new_mutation_smiles_list) == 0:
        new_ligands_list = []
    else:
        new_ligands_list = new_mutation_smiles_list

    loop_counter = 0

    number_of_processors = int(vars["parallelizer"].return_node())

    # initialize the smileclickclass
    a_smiles_click_chem_object = SmileClickClass.SmilesClickChem(
        rxn_library_variables, new_mutation_smiles_list, vars["filter_object_dict"]
    )

    while loop_counter < 2000 and len(new_ligands_list) < num_mutants_to_make:

        react_list = copy.deepcopy(ligands_list)

        while len(new_ligands_list) < num_mutants_to_make and len(react_list) > 0:

            a_smiles_click_chem_object.update_list_of_already_made_smiles(new_ligands_list)
            num_to_grab = num_mutants_to_make - len(new_ligands_list)
            num_to_make = num_to_grab

            # to minimize a big loop of running a single mutation at a time we
            # will make 1 new lig/processor. This will help to prevent wasting
            # reasources and time. 为了最大限度地减少一次运行单个突变的大循环，我们将制作 1 个新的 lig/processor。这将有助于防止浪费资源和时间。
            if num_to_make < number_of_processors:
                num_to_make = number_of_processors

            smile_pairs = [
                react_list.pop() for x in range(num_to_make) if len(react_list) > 0
            ]

            smile_inputs = [x[0] for x in smile_pairs]
            smile_names = [x[1] for x in smile_pairs]

            job_input = tuple(
                [tuple([smile, a_smiles_click_chem_object]) for smile in smile_inputs]
            )

            results = vars["parallelizer"].run(
                job_input, run_smiles_click_for_multithread
            )

            for index, i in enumerate(results):
                if i is not None:

                    # Get the new molecule's (aka the Child lig) Smile string 获取新分子（又名 Child lig）的 Smile 字符串
                    child_lig_smile = i[0]

                    # get the reaction id number 获取反应 ID 号
                    reaction_id_number = i[1]

                    # get the ID for the parent of a child mol and the
                    # complementary parent mol. comp mol could be None or a
                    # zinc database ID 获取子 mol 的父级和互补父级 mol 的 ID。comp mol 可以是 None 或zinc数据库 ID
                    parent_lig_id = smile_names[index]
                    zinc_id_comp_mol = i[2]

                    # Make a list of all smiles and smile_id's of all
                    # previously made smiles in this generation 列出这一代中所有微笑以及之前所有smile的 smile_id
                    list_of_already_made_smiles = []
                    list_of_already_made_id = []

                    # fill lists of all smiles and smile_id's of all
                    # previously made smiles in this generation 填充这一代中所有微笑和所有先前smile的 smile_id 列表
                    for x in new_ligands_list:
                        list_of_already_made_smiles.append(x[0])
                        list_of_already_made_id.append(x[1])

                    if child_lig_smile not in list_of_already_made_smiles:
                        # if the smiles string is unique to the list of 如果smile字符串对于本轮反应中先前smile字符串的列表是唯一的，
                        # previous smile strings in this round of reactions
                        # then we append it to the list of newly created
                        # ligands we append it with a unique ID, which also
                        # tracks the progress of the reactant 那么我们将其附加到新创建的配体列表中，我们将其附加一个唯一的ID，该ID还跟踪反应物的进度
                        is_name_unique = False
                        while is_name_unique is False:
                            # make unique ID with the 1st number being the
                            # parent_lig_id for the derived mol, Followed by
                            # Mutant, folowed by the generationnumber,
                            # followed by a unique. 
                            # 制作唯一 ID，第一个数字是派生摩尔的parent_lig_id，后面是突变体，后面是世代号，后面是唯一的。
                            # get the unique ID (last few diget ID of the
                            # parent mol
                            parent_lig_id = parent_lig_id.split(")")[-1]

                            random_id_num = random.randint(100, 1000000)
                            if zinc_id_comp_mol is None:
                                new_lig_id = "({})Gen_{}_Mutant_{}_{}".format(
                                    parent_lig_id,
                                    generation_num,
                                    reaction_id_number,
                                    random_id_num,
                                )
                            else:
                                new_lig_id = "({}+{})Gen_{}_Mutant_{}_{}".format(
                                    parent_lig_id,
                                    zinc_id_comp_mol,
                                    generation_num,
                                    reaction_id_number,
                                    random_id_num,
                                )

                            # check name is unique
                            if new_lig_id not in list_of_already_made_id:
                                is_name_unique = True

                        # make a temporary list containing the smiles string
                        # of the new product and the unique ID 制作一个临时列表，其中包含新产品的微笑字符串和唯一ID
                        ligand_info = [child_lig_smile, new_lig_id]

                        # append the new ligand smile and ID to the list of
                        # all newly made ligands 将新的配体smile和ID添加到所有新制作的配体的列表中
                        new_ligands_list.append(ligand_info)

        loop_counter = loop_counter + 1

    if len(new_ligands_list) < num_mutants_to_make:
        return None

    # once the number of mutants we need is generated return the list
    return new_ligands_list


def run_smiles_click_for_multithread(smile, a_smiles_click_chem_object):
    """
    This function takes a single smilestring and performs SmileClick on it.
    此函数采用单个smile字符串并对其执行 SmileClick。
    This is necessary for Multithreading as it is unable to execute
    multithread on a class function, but can thread a class run within a
    function.
    这对于多线程是必要的，因为它无法在类函数上执行多线程，但可以在函数内线程化运行的类。
    Inputs:
    :param str smile: a SMILES string

    Returns:
    :returns: str result_of_run: either a smile string of a child mol or None
        if the reactions failed
    """

    result_of_run = a_smiles_click_chem_object.run_smiles_click(smile)

    return result_of_run
