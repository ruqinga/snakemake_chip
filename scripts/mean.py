import pandas as pd
import argparse

def group_and_average(input_file, group_cols, avg_col, output_file):
    """
    按指定列分组计算某列的平均值，并输出为新文件。

    参数：
    - input_file: str, 输入的 txt 文件路径
    - group_cols: list, 用于分组的列索引（从 1 开始计数）
    - avg_col: int, 需要计算平均值的列索引（从 1 开始计数）
    - output_file: str, 输出的文件路径
    """
    # 读取文件，假定无列名
    print(f"Reading file: {input_file}")
    df = pd.read_csv(input_file, sep="\t", header=None)

    # 动态生成列名
    num_columns = df.shape[1]
    column_names = [f"col{i + 1}" for i in range(num_columns)]
    df.columns = column_names

    # 打印列名供检查
    print("Generated column names:")
    print(column_names)

    # 将列索引转换为列名
    group_col_names = [f"col{col}" for col in group_cols]
    avg_col_name = f"col{avg_col}"

    # 检查列索引是否超出范围
    for col in group_col_names + [avg_col_name]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found. Check column indices.")

    # 分组计算平均值
    print(f"Grouping by columns: {group_col_names} and calculating mean for '{avg_col_name}'")
    result = df.groupby(group_col_names)[avg_col_name].mean().reset_index()

    # 输出文件
    result.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    # 定义命令行参数
    parser = argparse.ArgumentParser(description="Group by specific columns and calculate the mean of another column.")
    parser.add_argument("input_file", help="Input txt file path")
    parser.add_argument("output_file", help="Output txt file path")
    parser.add_argument("-g", "--group", required=True,
                        help="Columns to group by (comma-separated indices starting from 1)")
    parser.add_argument("-v", "--value", required=True, type=int,
                        help="Column to calculate the average (index starting from 1)")

    # 解析参数
    args = parser.parse_args()

    # 处理分组列
    group_cols = [int(col) for col in args.group.split(",")]

    # 调用函数
    group_and_average(args.input_file, group_cols, args.value, args.output_file)
