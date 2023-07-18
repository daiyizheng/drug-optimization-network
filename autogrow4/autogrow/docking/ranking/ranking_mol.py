"""
This script runs the ranking and selection of ligands.
"""
import __future__

import os
import random

import rdkit
import rdkit.Chem as Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

import autogrow.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
import autogrow.docking.ranking.selecting.rank_selection as Rank_Sel
import autogrow.docking.ranking.selecting.roulette_selection as Roulette_Sel
import autogrow.docking.ranking.selecting.tournament_selection as Tournament_Sel


def create_seed_list(usable_list_of_smiles, num_seed_diversity,
                     num_seed_dock_fitness, selector_choice, tourn_size):
    """
    this function will take ausable_list_of_smiles which can be derived from
    either the previous generation or a source_compounds_file. Then it will
    select a set of smiles chosen by a weighted function to their
    docking-fitness (docking score) Then it will select a set of smiles chosen
    by a weighted function to their diversity-fitness (docking score) Then it
    will merge these sets of smiles into a single list
    该函数将采用 ausable_list_of_smiles，它可以从上一代或 source_compounds_file 派生。然后它将选择由加权函数选择的一组微笑到其对接适应度（对接分数） 然后它将选择由加权函数选择的一组微笑到其多样性适应度（对接分数） 然后它将合并这些集合将微笑放入一个列表中
    Using the merged list it will make a list of all the smiles in the
    merged-list (chosen_mol_list) with all the other information about each of
    the chosen mols from the usable_list_of_smiles
    使用合并列表，它将创建合并列表（chosen_mol_list）中所有微笑的列表，以及有关来自 usable_list_of_smiles 的每个所选摩尔的所有其他信息
    It will return this list with the complete information of each chosen mol
    (weighted_order_list)
    它将返回此列表，其中包含每个所选摩尔的完整信息（weighted_order_list）
    Inputs:
    :param list usable_list_of_smiles: a list with SMILES strings, names, and
        information about the smiles from either the previous generation or the
        source compound list
    :param int num_seed_diversity: the number of seed molecules which come
        from diversity selection
    :param int num_seed_dock_fitness: the number of seed molecules which come
        from eite selection by docking score
    :param int selector_choice: the choice of selector method. Choices are
        Roulette_Selector, Rank_Selector, or Tournament_Selector
    :param float tourn_size: percentage of the total pool of ligands to be
        tested in each tournament.

    Returns:
    :returns: list chosen_mol_full_data_list: a list of all the smiles in a
        weighted ranking ie ["CCCC"  "zinc123"   1    -0.1]
    """

    if selector_choice == "Roulette_Selector":

        print("Roulette_Selector")
        # Get seed molecules based on docking scores 根据对接分数获取种子分子
        docking_fitness_smiles_list = Roulette_Sel.spin_roulette_selector(
            usable_list_of_smiles, num_seed_dock_fitness, "docking"
        )

        # Get seed molecules based on diversity scores
        diversity_smile_list = Roulette_Sel.spin_roulette_selector(
            usable_list_of_smiles, num_seed_diversity, "diversity"
        )

    elif selector_choice == "Rank_Selector":
        print("Rank_Selector")
        # This assumes the most negative number is the best option which is
        # true for both This is true for both the diversity score and the
        # docking score. This may need to be adjusted for different scoring
        # functions. 这假设最大负数是最佳选择，这对于多样性分数和对接分数都是如此。这可能需要针对不同的评分函数进行调整。

        # Get seed molecules based on docking scores 根据对接分数获取种子分子
        docking_fitness_smiles_list = Rank_Sel.run_rank_selector(
            usable_list_of_smiles, num_seed_dock_fitness, -2, False
        )

        # Get seed molecules based on diversity scores 根据多样性分数获取种子分子
        diversity_smile_list = Rank_Sel.run_rank_selector(
            usable_list_of_smiles, num_seed_diversity, -1, False
        )

    elif selector_choice == "Tournament_Selector":
        print("Tournament_Selector")
        # This assumes the most negative number is the best option which is
        # true for both This is true for both the diversity score and the
        # docking score. This may need to be adjusted for different scoring
        # functions.

        # Get seed molecules based on docking scores
        docking_fitness_smiles_list = Tournament_Sel.run_Tournament_Selector(
            usable_list_of_smiles, num_seed_dock_fitness, tourn_size, -2, True
        )

        # Get seed molecules based on diversity scores
        diversity_smile_list = Tournament_Sel.run_Tournament_Selector(
            usable_list_of_smiles, num_seed_diversity, tourn_size, -1, True
        )

    else:
        print(selector_choice)
        raise Exception(
            "selector_choice value is not Roulette_Selector, Rank_Selector, nor Tournament_Selector"
        )

    chosen_mol_list = [x for x in docking_fitness_smiles_list]
    chosen_mol_list.extend(diversity_smile_list)

    if selector_choice in ["Rank_Selector", "Roulette_Selector"]:
        # Get all the information about the chosen molecules. chosen_mol_list 获取有关所选分子的所有信息。
        # is 1D list of all chosen ligands chosen_mol_full_data_list is a 1D 是所有选定配体的一维列表 selected_mol_full_data_list 是一个一维列表，
        # list with each item of the list having multiple pieces of
        # information such as the ligand name/id, the smiles string, the
        # diversity and docking score...  selected_mol_list 列表中的每个项目都具有多个信息，例如配体name/id、smiles字符串、多样性和对接分数...
        chosen_mol_full_data_list = get_chosen_mol_full_data_list(
            chosen_mol_list, usable_list_of_smiles
        )

    elif selector_choice == "Tournament_Selector":
        # Tournament_Selector returns an already full list of ligands so you
        # can skip the get_chosen_mol_full_data_list step
        chosen_mol_full_data_list = chosen_mol_list

    else:
        print(selector_choice)
        raise Exception(
            "selector_choice value is not Roulette_Selector, Rank_Selector, nor Tournament_Selector"
        )

    return chosen_mol_full_data_list


def get_chosen_mol_full_data_list(chosen_mol_list, usable_list_of_smiles):
    """
    This function will take a list of chosen molecules and a list of all the
    SMILES which could have been chosen and all of the information about those
    SMILES (ie. ligand name, SMILES string, docking score, diversity score...)
    该函数将获取所选分子的列表以及可以选择的所有 SMILES 的列表以及有关这些 SMILES 的所有信息（即配体名称、SMILES 字符串、对接分数、多样性分数...）
    It will iterated through the list of chosen mols (chosen_mol_list), get
    all the information from the usable_list_of_smiles Then it appends the
    corresponding item in usable_list_of_smiles to a new list
    weighted_order_list
    它将迭代所选摩尔列表（chosen_mol_list），从 usable_list_of_smiles 中获取所有信息，然后将 usable_list_of_smiles 中的相应项目附加到新列表weighted_order_list
    --- an issue to be aware of is that there may be redundancies in both
        chosen_mol_list and usable_list_of_smiles this causes a many-to-many
        problem so if manipulating this section you need to solve for
        one-to-many
    ---for this reason if this gets altered it will raise an
        AssertionError if the one-to-many is violated.

    It then shuffles the order of the list which to prevent biasing by the
    order of the ligands.

    It will return that list of the chosen molecules in a randomly shuffled
    order.

    Inputs:
    :param list chosen_mol_list: a list of chosen molecules
    :param list usable_list_of_smiles: List of all the possibly chosen ligs
        and all the of the info about it (ie. ligand name, SMILES string, docking
        score, diversity score...) ["CCCC"  "zinc123"   1    -0.1  -0.1]

    Returns:
    :returns: list weighted_order_list: a list of all the SMILES with all of
        the associated information in a random order
    """

    sorted_list = sorted(usable_list_of_smiles, key=lambda x: float(x[-2]))
    weighted_order_list = []
    for smile in chosen_mol_list:
        for smile_pair in sorted_list:
            if smile == smile_pair[0]:
                weighted_order_list.append(smile_pair)
                break

    if len(weighted_order_list) != len(chosen_mol_list):
        raise AssertionError(
            "weighted_order_list not the same length as the chosen_mol_list"
        )

    random.shuffle(weighted_order_list)

    return weighted_order_list


def get_usable_format(infile):
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.
    此代码为格式化为.smi文件的文件获取一个字符串。它打开文件并将组件读入可用列表中。
    The .smi must follow the following format for each line:.smi对于每一行都必须遵循以下格式：
        MANDATORY INFO 强制性信息
            part 1 is the SMILES string  第1部分是SMILES字符串
            part 2 is the SMILES name/ID 第2部分是SMILES名称/ID

        Optional info 可选信息
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population 第1部分（最后一条信息）是SMILES相对于其种群的多样性得分
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating 第2部分（倒数第二条信息）是用于评估的适应度度量
                - For default setting this is the Docking score 对于默认设置，这是对接分数
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2 如果添加一个唯一的评分函数，则对接分数应为-3，该评分函数应为-2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)
            任何其他信息必须在第2部分和第2部分之间（这允许在不中断代码其余部分的情况下扩展功能）
    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles


def convert_usable_list_to_lig_dict(usable_list_of_smiles):
    """
    This will convert a list created by get_usable_format() to a dictionary
    using the ligand smile+lig_id as the key. This makes for faster searching
    in larger Autogrow runs.

    Inputs:
    :param list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow

    Returns:
    :returns: list usable_dict_of_smiles: djct of all the ligand info with a
        key containing both the SMILE string and the unique lig ID
    """

    if type(usable_list_of_smiles) is not type([]):
        return None

    usable_dict_of_smiles = {}
    for item in usable_list_of_smiles:
        key = str(item[0]) + str(item[1])
        if key in usable_dict_of_smiles.keys():
            if usable_dict_of_smiles[key][-1] < item[-2]:
                continue
        usable_dict_of_smiles[key] = item
    return usable_dict_of_smiles


##### Called in the docking class ######
def score_and_append_diversity_scores(molecules_list):
    """
    This function will take list of molecules which makes up a population. It
    will then create a diversity score for each molecules:
    It creates the diversity score by determining the Morgan Fingerprint for
        each molecule in the population.
    It then compares the fingerprints for every molecule against every
        molecule in a pairwise manner.
        Based on the approach provided on
            http://www.rdkit.org/docs/GettingStartedInPython.html section: "Picking
            Diverse Molecules Using Fingerprints"
        It determines a score of similarity using the RDKit function
            DataStructs.DiceSimilarity
            -The higher the similarity the higher the similarity score
                -ie) if you compare two identical SMILES the similarity score
                    is 1.0. I.e., if you compare 4 identical SMILES the
                    similarity score for each is 4.0.
                -ie) if you compare two completely different SMILES, the score
                    is 0.0

        It sums the similarity score for each pairwise comparison.
            -ie) if there are 15 ligands the max score is 15 the minimum is 0.0
                    with 15.0 if all ligands are identical

        It then appends the diversity score to the molecule list which it
        returns.

        It can raise an AssertionError if there are ligs which fail to
            sanitize or deprotanate.
                -this prevents future errors from occuring in later steps and
                    makes this funciton usable for multiple codes
        It will remove any Nones from the input list

    Inputs:
    :param list molecules_list: list of all molecules in the populations with
    the respective info

    Returns:
    :returns: list molecules_list: list of all molecules in the populations
        with the respective info and append diversity score
    """

    mol_list = []

    for pair in molecules_list:
        if pair is not None:
            smile = pair[0]
            # name = pair[1]
            try:
                mol = Chem.MolFromSmiles(smile, sanitize=False)
            except:
                mol = None

            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                    def score_and_append_diversity_scores"
                )

            mol = MOH.check_sanitization(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def score_and_append_diversity_scores"
                )

            mol = MOH.try_deprotanation(mol)
            if mol is None:
                raise AssertionError(
                    "mol in list failed to sanitize. Issue in Ranking.py \
                                        def score_and_append_diversity_scores"
                )

            temp = [x for x in pair]
            temp.append(mol)
            if temp[-1] is None:
                print(temp)
                print("None in temp list, skip this one")
                continue
            if temp[-1] is not None:
                mol_list.append(temp)
        else:
            print("noneitem in molecules_list in score_and_append_diversity_scores")

    fps_list = []
    for molecule in mol_list:
        fp = GetMorganFingerprint(molecule[-1], 10, useFeatures=True) # 计算指纹
        temp = [x for x in molecule]
        temp.append(fp)
        fps_list.append(temp)

    fps_list_w_div_score = []
    for i in range(0, len(fps_list)):
        diversity_score = 0
        for j in range(0, len(fps_list)):
            if i != j:
                # if DiceSimilarity=1.0 its a perfect match, the smaller the
                # number the more diverse it is. The sum of all of these gives
                # the distance from the normal. The smaller the number means
                # the more distant 如果 DiceSimilarity=1.0 则为完美匹配，数字越小则越多样化。所有这些的总和给出了与法线的距离。数字越小代表距离越远
                diversity_score = diversity_score + DataStructs.DiceSimilarity(
                    fps_list[i][-1], fps_list[j][-1]
                )
        temp = [x for x in fps_list[i]]
        temp.append(str(diversity_score))
        fps_list_w_div_score.append(temp)

    # take the diversity score and append to the last column in the original
    # list 获取多样性分数并附加到原始列表中的最后一列

    for i in range(0, len(molecules_list)):
        if molecules_list[i][0] == fps_list_w_div_score[i][0]:

            molecules_list[i].append(fps_list_w_div_score[i][-1])

    return molecules_list
