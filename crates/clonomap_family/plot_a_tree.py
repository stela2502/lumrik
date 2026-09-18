from ete3 import Tree, TreeStyle

t = Tree("/home/med-sal/sens05_shared/jyuan/no_backup/GM/Stefans_analysis/ChangeO_Db_2025/DefineClones_2025/A_2025_ProductiveCloneDfined_clone_123_len_385_input_tree.newick")

ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 200    # spacing between nodes

t.show(tree_style=ts)
