#!/bin/bash
# $1: input TSV file
# $2: output Phylip distance matrix file

# 获取样本数量（不包括标题行）
n_samples=$(tail -n +2 "$1" | wc -l)
echo $n_samples > "$2"

# 处理每一行：第一列作为样本名称，其余列作为距离值
tail -n +2 "$1" | awk '{
    # 输出样本名称
    printf "%s", $1
    # 输出距离值（从第二列开始）
    for (i=2; i<=NF; i++) {
        printf " %.15g", $i
    }
    print ""
}' >> "$2"
