import sys
from Bio import Phylo
from io import StringIO
from collections import defaultdict

def get_leaf_names(clade):
    """获取一个分支的所有叶子节点名称"""
    return frozenset(leaf.name for leaf in clade.get_terminals())

def build_support_dict(tree):
    """构建支持率字典，键为叶子节点集合，值为支持率"""
    support_dict = {}
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.confidence is not None:
            leaf_names = get_leaf_names(clade)
            support_dict[leaf_names] = str(clade.confidence)
        elif not clade.is_terminal() and clade.name is not None and clade.name.isdigit():
            leaf_names = get_leaf_names(clade)
            support_dict[leaf_names] = clade.name
    return support_dict

def add_support_to_tree(tree, support_dict):
    """将支持率添加到树中"""
    for clade in tree.find_clades():
        if not clade.is_terminal():
            leaf_names = get_leaf_names(clade)
            support = support_dict.get(leaf_names, "0")
            if clade.name is None:
                clade.name = support
            else:
                clade.name = f"{support}_{clade.name}"

def main():
    if len(sys.argv) != 4:
        print("Usage: python merge.py input_bootstrap input_length output")
        sys.exit(1)
    
    bootstrap_file = sys.argv[1]
    length_file = sys.argv[2]
    output_file = sys.argv[3]
    
    try:
        # 读取支持率树
        bootstrap_tree = next(Phylo.parse(bootstrap_file, 'newick'))
    except Exception as e:
        print(f"Error reading bootstrap tree: {e}")
        sys.exit(1)
    
    try:
        # 读取分支长度树
        length_tree = next(Phylo.parse(length_file, 'newick'))
    except Exception as e:
        print(f"Error reading length tree: {e}")
        sys.exit(1)
    
    # 构建支持率字典
    support_dict = build_support_dict(bootstrap_tree)
    
    # 将支持率添加到分支长度树
    add_support_to_tree(length_tree, support_dict)
    
    # 写入输出文件
    Phylo.write(length_tree, output_file, 'newick')
    print(f"Merged tree written to {output_file}")

if __name__ == "__main__":
    main()